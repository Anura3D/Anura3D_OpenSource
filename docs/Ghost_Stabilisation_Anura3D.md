# Ghost stabilisation of the explicit MPM in Anura3D

**Status:** implemented for 2D TRI3, explicit dynamic, dry single phase (`src/GhostStabilisation.FOR`).
Section 8 lists what differs from the original design. Parameter values still need calibrating.
**Scope:** 2D plane strain and 2D axisymmetric on TRI3 meshes; 3D TET4 noted where it differs.
Explicit dynamic MPM-MIXED, dry single-phase material (`NFORMULATION = 1`, one phase), with or
without Bardenhagen contact and the rigid body.
**Reference formulation:** W.M. Coombs (2023), *Ghost stabilisation of the Material Point Method for
stable quasi-static and dynamic analysis of large deformation problems*, IJNME 124(21):4841–4875,
doi:10.1002/nme.7332 (preprint arXiv:2209.10955). The method itself is Burman's ghost penalty.
The mass-matrix ghost term in explicit time stepping follows the idea of Stoter et al. (2023), CMAME
(arXiv:2302.07019).

> **Verification flag.** The operators below follow the structure of Coombs (2023). Check the
> penalty scalings and recommended parameter values against §3–§5 of that paper before coding.
> Values marked *[to check]* are proposals, not quotations.

---

## 1. The problem in Anura3D terms

### 1.1 Where the instability comes from

In the explicit solver, the soil entity's nodal quantities are mapped from material points (MPs):

$$
m_a=\sum_p m_p\,N_a(\mathbf x_p),\qquad
\mathbf f^{\text{int}}_a=\sum_p V_p\,\boldsymbol\sigma_p\,\nabla N_a(\mathbf x_p),\qquad
\mathbf a_a=\frac{\mathbf f^{\text{ext}}_a-\mathbf f^{\text{int}}_a}{m_a}.
$$

Take a node $a$ that is reached only by MPs sitting close to the opposite side of their element.
Then $N_a(\mathbf x_p)\to0$, so $m_a\to0$. But $\nabla N_a$ is constant over a linear element and of
order $1/h$, so $\mathbf f^{\text{int}}_a$ does **not** go to zero. The resulting $\mathbf a_a$ is
large and spurious. In an implicit code the same effect shows up as the near-singular stiffness that
Coombs describes. In Anura's explicit code it shows up as:

* spurious accelerations and velocity spikes at the edge nodes of the material body;
* stress oscillations in MPs near free surfaces and behind the impactor;
* MPs being thrown out ("outside the mesh") or elements inverting;
* a stable time step that depends on how much of the element is filled. `CalculateCriticalTimeStep`
  (`MPMDYNConvPhase.FOR:2882`) uses only `ElementLMin / c` and cannot see this.

### 1.2 Where it can occur

MPM-MIXED already integrates fully filled elements with Gauss points, which is stable. MP-based
integration is used only where `IsParticleIntegration(IEl)` is `.true.`, set by
`CheckFillingOfElements` (`MPMConvPhase.FOR:952`):

* the element holds a boundary MP **and** $\sum_{p\in K}V_p \le$ `RequiredDegreeOfFilling`$\cdot|K|$; or
* the element holds MPs of more than one material.

Elements with very small, uneven fill are the source of the problem, so the stabilisation is
confined to these elements and their neighbours.

---

## 2. Notation

| Symbol | Meaning | Anura3D variable |
| --- | --- | --- |
| $K$, $\lvert K\rvert$ | element, its area (2D) / volume (3D) | `IEl`, `ElementSpace(IEl)` |
| $F$, $\lvert F\rvert$ | interior face (edge in 2D), its length / area | side `ISide` of `IEl` |
| $K^+,K^-$ | the two elements sharing $F$ | `IEl`, `ElementAdjacencies(IEl,ISide)` |
| $\mathbf n_F$ | unit normal of $F$, from $K^+$ to $K^-$ | computed from side nodes |
| $h_F$ | face length scale | `min(ElementLMin(K+),ElementLMin(K-))` |
| $N_a$, $\nabla N_a$ | linear shape function, constant gradient per element | from `FormB3` (`BuildBJacDet.FOR:40`) |
| $\phi_K$ | fill fraction, $\sum_{p\in K}V_p/\lvert K\rvert$ | `IntegrationWeight`, `ElementSpace` |
| $m_a$ | lumped soil mass at node $a$ | `LumpedMassDry(iDof, SOFT_ENTITY)` |
| $\mathbf v_a^n$, $\Delta t$ | nodal velocity from mapping, time increment | `TotalVelocitySoil`, `CalParams%TimeIncrement` |
| $\rho$, $c$, $M^*=\rho c^2$ | density, P-wave speed, P-wave modulus | `GetWaveSpeedForTimeStepSize` |
| $\varrho_F$ | integration weight: $1$ (plane strain) or $2\pi r_F$ (axisymmetric) | radius of face midpoint |
| $\gamma_K,\gamma_M$ | ghost stiffness and mass penalty parameters | new CPS input |
| $[\![\cdot]\!]$ | jump across $F$: $(\cdot)\vert_{K^+}-(\cdot)\vert_{K^-}$ | |

---

## 3. Ghost faces

The ghost face set $\mathcal F_G$ contains every interior face $F=\bar K^+\cap\bar K^-$ such that:

1. both $K^+$ and $K^-$ are active (`IsActiveElement`) and contain soil-entity MPs;
2. at least one of them is **cut**: $\phi_K<\phi_G$ ($\phi_G$ is a new input, default 0.6).
   `IsParticleIntegration` is not used, because it follows the 0.9 degree of filling. On the cone
   model the fill of interior elements fluctuates between about 0.7 and 1.05 as MPs migrate. A
   threshold of 0.9 flagged about half of the interior; 0.6 flags about 1 % of elements, the actual
   edge elements;
3. neither element contains rigid-body MPs (`EntityIDArray == CalParams%RigidBody%RigidEntity`);
4. neither element contains only gas MPs (when gasification is active; see `Gasification_Anura3D.md`).

Mesh faces on the outer mesh boundary are never ghost faces, because `ElementAdjacencies` returns
`0` there. The list is rebuilt every time step, because the set of cut elements changes as MPs move.
The neighbour search already exists: `DetermineElementSideAdjacencies` (`ElemConnections.FOR:824`)
fills `ElementAdjacencies(IEl, ISide)` for TRI3 (sides 1-2, 3-1, 2-3) and for tetrahedra.

Each face is stored once (for example only when `IEl < ElementAdjacencies(IEl,ISide)`), with its two
elements, its 4 (TRI3) or 5 (TET4) patch nodes, $\lvert F\rvert$, $\mathbf n_F$, $h_F$ and $\varrho_F$.

---

## 4. Ghost operators for linear elements

For polynomial order $p=1$ only the jump in the first normal derivative is penalised (Burman's ghost
penalty with $j=1$). Because $\nabla N_a$ is constant per element, the jump is exact and constant
over the face:

$$
[\![\partial_n N_a]\!]_F=\big(\nabla N_a\vert_{K^+}-\nabla N_a\vert_{K^-}\big)\cdot\mathbf n_F .
$$

It is non-zero only for the nodes of $K^+\cup K^-$. A node that belongs to only one of the two
elements has $\nabla N_a=0$ in the other.

### 4.1 Ghost stiffness

$$
K^g_{ab}=\gamma_K\,M^*_F\sum_{F\in\mathcal F_G} h_F\,\varrho_F\,\lvert F\rvert\;
[\![\partial_n N_a]\!]_F\,[\![\partial_n N_b]\!]_F\;\mathbf I_{d\times d},
\qquad M^*_F=\max(\rho c^2)_{K^+,K^-}.
$$

**Dimension check:** $[M^*]\,[L]\,[L^{d-1}]\,[L^{-2}]=[M^*][L^{d-2}]$, which matches
$K=\int B^{\mathsf T}DB\,\mathrm d\Omega$.

### 4.2 Ghost mass

$$
M^g_{ab}=\gamma_M\,\rho_F\sum_{F\in\mathcal F_G} h_F^{3}\,\varrho_F\,\lvert F\rvert\;
[\![\partial_n N_a]\!]_F\,[\![\partial_n N_b]\!]_F .
$$

**Dimension check:** $[\rho][L^3][L^{d-1}][L^{-2}]=[\rho][L^d]$, which matches
$M=\int\rho NN\,\mathrm d\Omega$. Here $\rho_F$ is the current mean MP density of the two elements.

### 4.3 Properties used below

1. **Affine fields are not penalised.** For any field $\mathbf w$ that is linear across $K^+\cup K^-$,
   $[\![\partial_n\mathbf w]\!]=0$. So $K^g$ and $M^g$ do not act on rigid translation, rigid
   rotation or homogeneous strain. They only damp kinks in the field across the ghost faces.
2. **Zero row sums.** $\sum_a N_a\equiv1$ gives $\sum_a\nabla N_a=0$ in each element, so
   $\sum_a[\![\partial_n N_a]\!]_F=0$. Therefore $\sum_a K^g_{ab}=\sum_a M^g_{ab}=0$.
3. **Symmetric positive semi-definite.** Both operators are sums of rank-one terms $\mathbf g\mathbf g^{\mathsf T}$
   with positive weights.

---

## 5. Adaptation to the explicit solver

### 5.1 Ghost force (stiffness term)

Anura3D is an updated-Lagrangian MPM: the grid is reset every step, so there is no total nodal
displacement for $K^g$ to act on. The ghost stiffness therefore acts on the **step increment**,
evaluated with the nodal velocities mapped at the start of the step:

$$
\Delta\mathbf u_b=\Delta t\,\mathbf v^n_b,\qquad
\mathbf f^g_a=\sum_b K^g_{ab}\,\Delta\mathbf u_b=\Big(\Delta t\,K^g\Big)_{ab}\mathbf v^n_b .
$$

This is a **ghost viscosity** $C^g=\Delta t\,K^g$. It is added to the internal force of the soil entity:

$$
\mathbf f^{\text{int}}_a\leftarrow\mathbf f^{\text{int}}_a+\mathbf f^g_a .
$$

* **Momentum:** property 2 gives $\sum_a\mathbf f^g_a=0$, so total linear momentum is conserved
  **exactly**. The ghost force only redistributes momentum inside each face patch.
* **Energy:** $\mathbf v^{\mathsf T}C^g\mathbf v\ge0$, so the term is dissipative and only removes
  energy from kinked (non-affine) velocity modes.
* **Contact:** it acts on `TotalVelocitySoil(:,SOFT_ENTITY)` only. The rigid entity is excluded
  (rule 3 in §3), and no face crosses the soil–body interface, so the Bardenhagen correction is untouched.

### 5.2 Ghost mass

Row-sum lumping, which Anura uses for $M_L$, **cannot** be used for $M^g$: by property 2 every row
sum is zero, so a lumped $M^g$ vanishes identically. Two options are provided.

**Option A: Jacobi solve of the consistent ghost system (recommended).** Solve

$$
(M_L+M^g)\,\mathbf a=\mathbf f,\qquad
D=M_L+\operatorname{diag}(M^g),\qquad
\mathbf a^{k+1}=D^{-1}\Big(\mathbf f-\big(M^g-\operatorname{diag}(M^g)\big)\mathbf a^{k}\Big),
\quad \mathbf a^0=D^{-1}\mathbf f ,
$$

with $k_{\max}=2$–$3$ sweeps over the ghost-patch nodes only. All other nodes keep
$\mathbf a=\mathbf f/m_a$. At convergence $\sum_a m_a\mathbf a_a=\sum_a\mathbf f_a$, because the ghost
rows sum to zero, so momentum is conserved. With a finite number of sweeps the residual momentum
error is logged (§7).

**Option B: diagonal ghost mass (mass scaling of small nodes).**
$\tilde m_a=m_a+M^g_{aa}$ is used in the acceleration only. Nodal velocities are still computed from
mapped momentum and true mass in the convective phase. This is cheaper but **not** momentum
conserving: the per-step error is $\sum_a M^g_{aa}\mathbf a_a$. It is useful for a first
implementation and as a fallback. A cap $M^g_{aa}\le\beta\,\bar m$ (with $\bar m$ the mean nodal
mass of the patch) bounds the error.

**Where:** `CalculateIncrementalNodalAcceleration` (`LagrangianPhase.FOR:1447`). For the soil entity,
replace `AccelerationSoil = RateofMomentum / LumpedMassDry` on ghost-patch DOFs by option A or B. The
system acceleration `AccelerationSys` used by the contact algorithm must use the same corrected
masses on those DOFs (`LumpedMassSys` plus the ghost diagonal), so that contact sees a consistent
soil–system pair.

### 5.3 Critical time step

The ghost mass raises small $m_a$, which allows a larger time step. The ghost stiffness adds
stiffness and damping, which require a smaller one. Add a nodal check for ghost-patch nodes in
`CalculateCriticalTimeStep` (`MPMDYNConvPhase.FOR:2882`):

* **Stiffness (Gershgorin) bound:**
  $\Delta t_a\le C\sqrt{\tilde m_a\big/\textstyle\sum_b\lvert K^{\text{el}}_{ab}+K^g_{ab}\rvert}$.
  Here $K^{\text{el}}_{ab}\approx M^*\,\lvert K\rvert\,\nabla N_a\!\cdot\!\nabla N_b$ is summed over
  the patch, $\tilde m_a=m_a+M^g_{aa}$, and $C$ is `CalParams%CourantNumber`.
* **Viscous bound** for the explicit ghost viscosity: $\Delta t\,c^g_{aa}/\tilde m_a\le 2$, where
  $c^g_{aa}=\Delta t\,K^g_{aa}$. This gives
  $\Delta t^2\,K^g_{aa}/\tilde m_a\le2$, a direct upper limit on $\gamma_K$ for a given time step,
  which is reported in the log.

The global step is the minimum of the existing element step and these nodal steps.

### 5.4 Parameters

| Parameter | Meaning | Proposed start | Range to study |
| --- | --- | --- | --- |
| $\gamma_K$ | ghost stiffness penalty | $0.05$ *[to check]* | $10^{-3}$ – $1$ |
| $\gamma_M$ | ghost mass penalty | $0.05$ *[to check]* | $10^{-3}$ – $1$ |
| $\phi_G$ | fill below which an element is "cut" | 0.6 | 0.4 – 0.7 |
| $k_{\max}$ | Jacobi sweeps (option A) | 2 | 1 – 5 |

Coombs (2023) reports results that are insensitive to the penalty over several orders of magnitude
in the implicit setting. In explicit form, $\gamma_K$ is also bounded by the viscous limit in §5.3.

---

## 6. Implementation plan

### 6.1 New module `src/GhostStabilisation.FOR` (`ModGhostStabilisation`)

| Routine | Purpose |
| --- | --- |
| `InitialiseGhost()` / `DestroyGhost()` | allocate and free face arrays (size $\le$ `NEl*ELEMENTSIDES/2`) |
| `BuildGhostFaces()` | apply the rules of §3, store patch nodes, $\lvert F\rvert$, $\mathbf n_F$, $h_F$, $\varrho_F$ and the jump vectors $[\![\partial_n N]\!]_F$ |
| `GhostForce(Vel, Fg)` | $\mathbf f^g=\Delta t\,K^g\mathbf v$, assembled face by face (never forms a global matrix) |
| `GhostMassDiag(Mg)` | $M^g_{aa}$ for option B and for the Jacobi preconditioner |
| `GhostJacobiSolve(F, M, A)` | option A, applied only to ghost-patch DOFs |
| `GhostTimeStep(dt)` | nodal bounds of §5.3 |
| `GhostDiagnostics()` | number of faces, $\max M^g_{aa}/m_a$, momentum error |

### 6.2 Hooks in the existing code

| Step | File / location | Change |
| --- | --- | --- |
| Rebuild faces | `MPMDYNConvPhase.FOR` after `SetActiveElement` / `CheckFillingOfElements` (≈ l.490–493) | `call BuildGhostFaces()` |
| Ghost force | `LagrangianPhase.FOR` after `call MPMDYNBTSig(...)` (l.1036) | `IntLoad(:,SOFT_ENTITY) += Fg` |
| Ghost mass | `CalculateIncrementalNodalAcceleration` (`LagrangianPhase.FOR:1447`) | option A or B on ghost DOFs, soil and system |
| Time step | `CalculateCriticalTimeStep` (`MPMDYNConvPhase.FOR:2882`) | `min` with `GhostTimeStep` |
| Output | `WriteVTKOutput.FOR` | element scalar `ghost_element`, nodal `ghost_mass_ratio` |

Everything is guarded by `if (.not.CalParams%Ghost%Apply) return`, so results are bit-identical
when the feature is off.

### 6.3 Inputs and problem type

**Kernel.** New derived type in `CalParams` (`ReadCalculationData.FOR`):

```fortran
type GhostType
  logical :: Apply = .false.
  real(REAL_TYPE) :: GammaK = 0.05, GammaM = 0.05, FillThreshold = 0.6
  integer(INTEGER_TYPE) :: MassSolver = 1   ! 0 = diagonal (B), 1 = Jacobi (A)
  integer(INTEGER_TYPE) :: NJacobi = 2
end type
```

* Defaults are set in `InitialiseCalculationParameters` with `Apply = .false.`, so old CPS files run
  unchanged.
* Reader branch next to `$$APIC_FLIP_BLEND_FACTOR` (≈ l.1826) for the key `$$GHOST_STABILISATION`, with
  `Assert` checks: $\gamma\ge0$, $0<\phi_G\le1$, `MassSolver` ∈ {0,1}, `NJacobi` ≥ 1.
* The kernel's CPS writer (`WriteCalculationParameters`, l.2887; `WriteReal` block ≈ l.3140) re-emits
  the key, so continuation files `CPS_002…` keep it. Add `$$APIC_FLIP_BLEND_FACTOR` there too; today
  it is missing from continuation files.
* Echo the active values to the log at stage start.

**CPS format** (one key, one line, as with `$$MASS_SCALING`):

```
$$GHOST_STABILISATION
1  0.05  0.05  0.60  1  2        ! flag  gammaK  gammaM  phiG  massSolver  nJacobi
```

**Problem type** (`src/GiD_Problemtype/Anura3D_2025.gid`). The new container is per stage, so
stabilisation can be off during gravity initialisation and on during impact:

```
Calculation Data
└─ STABILISATION AND GASIFICATION
   └─ GHOST STABILISATION   [do not apply ghost stabilisation | apply ghost stabilisation]
       ├─ ghost preset              [conservative | default | aggressive | user]
       ├─ ghost stiffness penalty [-]   (shown for "user")
       ├─ ghost mass penalty [-]        (shown for "user")
       ├─ ghost-face fill threshold [-]
       ├─ ghost mass solver         [diagonal | Jacobi]
       └─ Jacobi iterations [-]
```

`xml/CalculationData.xml`, copying the `MASS_SCALING` pattern (l.90–94):

```xml
<container n="STABILISATION" pn="STABILISATION AND GASIFICATION" icon="calculate.png">
  <value n="GHOST_STABILISATION" pn="GHOST STABILISATION" icon="calculate.png"
         values="do not apply ghost stabilisation,apply ghost stabilisation"
         v="do not apply ghost stabilisation" actualize_tree="1">
    <dependencies node="../value[@n='ghost_preset']" att1="state" v1="normal" value="apply ghost stabilisation"/>
    <dependencies node="../value[@n='ghost_preset']" att1="state" v1="hidden" value="do not apply ghost stabilisation"/>
    <!-- same pair for ghost_fill_threshold, ghost_mass_solver, ghost_jacobi_iterations -->
  </value>
  <value n="ghost_preset" pn="ghost preset" values="conservative,default,aggressive,user" v="default"
         state="hidden" actualize_tree="1"
         help="conservative: gK=gM=0.01; default: 0.05; aggressive: 0.2; user: enter values">
    <dependencies node="../value[@n='ghost_gamma_K']" att1="state" v1="normal" value="user"/>
    <dependencies node="../value[@n='ghost_gamma_K']" att1="state" v1="hidden" value="conservative"/>
    <!-- ... hidden for default/aggressive; same for ghost_gamma_M -->
  </value>
  <value n="ghost_gamma_K" pn="ghost stiffness penalty [-]" v="0.05" state="hidden"
         help="Penalty on gradient jumps across faces of cut elements. Range 1e-3 to 1."/>
  <value n="ghost_gamma_M" pn="ghost mass penalty [-]" v="0.05" state="hidden"
         help="Ghost mass added to small-mass nodes. Range 1e-3 to 1."/>
  <value n="ghost_fill_threshold" pn="ghost-face fill threshold [-]" v="0.9" state="hidden"
         help="Elements filled below this fraction are treated as cut. 0 to 1."/>
  <value n="ghost_mass_solver" pn="ghost mass solver" values="diagonal,Jacobi" v="Jacobi" state="hidden"/>
  <value n="ghost_jacobi_iterations" pn="Jacobi iterations [-]" v="2" state="hidden"/>
  <!-- GASIFICATION block: see Gasification_Anura3D.md -->
</container>
```

`scripts/createCPS.tcl`, placed next to `# APIC/FLIP BLEND FACTOR` (l.633) and modelled on
`# MASS SCALING` (l.240–250). Presets are resolved here, so the kernel only receives numbers:

```tcl
# GHOST STABILISATION
GiD_WriteCalculationFile puts {$$GHOST_STABILISATION}
set base {container[@n="Calculation_Data"]/container[@n="STABILISATION"]}
set on     [$stageNode selectNodes "string($base/value\[@n=\"GHOST_STABILISATION\"\]/@v)"]
set preset [$stageNode selectNodes "string($base/value\[@n=\"ghost_preset\"\]/@v)"]
switch $preset {
  conservative { set gK 0.01; set gM 0.01 }
  default      { set gK 0.05; set gM 0.05 }
  aggressive   { set gK 0.2;  set gM 0.2  }
  user { set gK [$stageNode selectNodes "string($base/value\[@n=\"ghost_gamma_K\"\]/@v)"]
         set gM [$stageNode selectNodes "string($base/value\[@n=\"ghost_gamma_M\"\]/@v)"] }
}
set phiG [$stageNode selectNodes "string($base/value\[@n=\"ghost_fill_threshold\"\]/@v)"]
set solv [expr {[$stageNode selectNodes "string($base/value\[@n=\"ghost_mass_solver\"\]/@v)"] eq "Jacobi" ? 1 : 0}]
set nJ   [$stageNode selectNodes "string($base/value\[@n=\"ghost_jacobi_iterations\"\]/@v)"]
set flag [expr {$on eq "apply ghost stabilisation" ? 1 : 0}]
GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s" $flag $gK $gM $phiG $solv $nJ]
```

---

## 7. Diagnostics and verification

**Implemented diagnostics.** At the first time step of each load step, the log file gets one line:
`Ghost: faces = …, max Mg/m = …, min m/mean m = …, ghost dt = …, momentum error = …`.
When the feature is on, the MP scalar `ghost_element` (1 = the MP's element touches a ghost face) is
added to the `MPScalar` VTK files.

**Per-step log line** (every `OUTPUT_CURVES_INTERVAL`):
number of ghost faces, $\max_a M^g_{aa}/m_a$, $\min_a m_a/\bar m$, controlling $\Delta t$ (element,
ghost stiffness or ghost viscous), and the momentum error
$\lVert\sum_a(m_a\mathbf a_a-\mathbf f_a)\rVert/\lVert\sum_a\mathbf f_a\rVert$.

**Benchmarks:**

1. **Elastic column under self-weight, offset on the grid** (as in Coombs 2023). Run the same column
   shifted by 0, ¼ and ½ element. Without ghost, the vertical stress profile oscillates and depends on
   the offset; with ghost it should match the analytical solution for every offset.
2. **Vibrating elastic bar, partially filled end element.** Check total energy (no growth), momentum
   error below $10^{-10}$ for option A, and a stable time step that does not depend on the fill fraction.
3. **Feature off.** Results must be bit-identical to the current build.
4. **`Free_fall_wall_BH` cone model.** Fewer velocity spikes at the free surface, no "outside the
   mesh" errors, and a cone reaction force within a few percent of the unstabilised run where that
   run is stable.
5. **Parameter sweep** over $\gamma_K,\gamma_M\in\{10^{-3},10^{-2},10^{-1},1\}$ to fix the presets.

---

## 8. As implemented (differences from sections 3–6)

| Topic | Implementation |
| --- | --- |
| Element types | 2D TRI3 only (plane strain and axisymmetric), `NFORMULATION = 1`, one phase. In any other case a log message is written and the feature switches itself off. |
| Ghost faces | Same rules as §3. Faces whose patch nodes include a contact interface node, or a node with mass of another entity, are also skipped. This keeps the rigid-body reaction force and the Bardenhagen correction unchanged. |
| Gradients | Computed analytically from the TRI3 node coordinates (`NodalCoordinatesUpd`), not from `FormB3`. |
| Axisymmetric weight | $\varrho_F = r_F$ (per radian), consistent with Anura's $r$-scaled MP weights. The fill fraction uses $V_p/r_p$. |
| Ghost force | Subtracted from `RateofMomentum(:,SOFT_ENTITY)` after damping, bulk viscosity and absorbing boundary, divided by `ScalingMassFactor` like every force (`LagrangianPhase.FOR`). |
| Ghost mass, Jacobi | Preconditioned Richardson iteration $\mathbf a\leftarrow\mathbf a+(\mathbf R-A\mathbf a)/P$ with $P_i = D_i+\sum_{j\ne i}\lvert M^g_{ij}\rvert$. $2P-A$ is diagonally dominant, so this converges for any mass distribution. Plain Jacobi can diverge. |
| **Momentum** | The consistent ghost system has affine modes with almost no mass and no penalty, so a few sweeps do not converge them. A per-cluster uniform acceleration correction is applied after the solve, $\mathbf a_i \mathrel{+}= \big(\sum_c R-\sum_c m\,a\big)/\sum_c m$, over the free DOFs of each connected ghost cluster $c$. A uniform field is affine, so it adds no ghost jump, and it restores $\sum m\,a=\sum R$ **exactly**, for both the diagonal and the Jacobi option. A check on a single face with masses of $10^{-8}$: momentum error $2\times10^{-16}$, max $\lvert a\rvert=0.88$, against $1.9\times10^7$ without ghost mass. |
| Time step | $\Delta t\le\sqrt{2\,s\,(m_a+M^g_{aa})/\sum_b\lvert K^g_{ab}\rvert}$ with mass scaling $s$, before the Courant factor (`CalculateCriticalTimeStep`). Keep $\gamma_M\gtrsim\gamma_K$, so this bound is about $\sqrt2\,h/c$ and not limiting. |
| CPS | `$$GHOST_STABILISATION` → `flag gammaK gammaM phiG massSolver nJacobi`. It is read in `ReadCalculationData.FOR` and re-written to continuation CPS files, together with `$$APIC_FLIP_BLEND_FACTOR`, which was missing before. |
| GUI | Container `STABILISATION` in `xml/CalculationData.xml` and a writer block in `scripts/createCPS.tcl`. Older models without the container write the key as off. |
