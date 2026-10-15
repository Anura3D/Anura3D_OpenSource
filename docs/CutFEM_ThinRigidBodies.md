# CutFEM treatment of elements cut by thin rigid bars in Anura3D

**Status:** design reference for review. No code has been written against this document yet.
**Scope:** 2D plane strain and 2D axisymmetric, TRI3 background mesh, explicit dynamic MPM.
**Target application:** cone penetration, where soil must be able to separate from the shaft.

---

## 1. Scope, notation and conventions

### 1.1 The problem

A thin rigid body in Anura3D is a chain of 2-node line ("bar") elements that moves through a
fixed background mesh of 3-node triangles. The bar chain is a *co-dimension-one* interface: it
has no volume, and it does not conform to the mesh. Elements it passes through are **cut**.

Today the background mesh carries a **single velocity field**, so all material points in a cut
element share the same nodal velocities and cannot separate. The purpose of this document is
to define a CutFEM-style discretisation in which they can.

### 1.2 Notation

| Symbol | Meaning |
| --- | --- |
| $\Omega$ | soil domain, split by the interface into $\Omega^+$ and $\Omega^-$ |
| $\Gamma$ | the interface — the bar chain |
| $\phi(\mathbf{x})$ | signed level set; $\phi>0$ in $\Omega^+$, $\phi<0$ in $\Omega^-$ |
| $\mathbf{n}$ | unit interface normal, oriented from $\Omega^-$ to $\Omega^+$ |
| $d$ | penetration, $d>0$ when a point has crossed into the body |
| $\mathbf{v}$, $\mathbf{a}$, $\boldsymbol\sigma$ | velocity, acceleration, Cauchy stress |
| $\mathbf{v}_R$ | rigid-body velocity field, $\mathbf{v}_R(\mathbf{x}) = \mathbf{v}_c + \omega \, \mathbf{e}_z \times (\mathbf{x}-\mathbf{x}_c)$ |
| $N_a$, $\mathbf{B}_a = \nabla N_a$ | shape function and its gradient at node $a$ |
| $h$ | element size; in code, `ElementLMin` (minimum altitude) |
| $K$, $F$ | an element and a face (edge) of the mesh |
| $[\![ \cdot ]\!]$ | jump across a face |
| subscript $p$, $a$ | material point, mesh node |

### 1.3 Plane strain versus axisymmetric

Every integral in this document carries a weight $\varrho$:

$$
\varrho = \begin{cases}
1 & \text{plane strain (per unit out-of-plane depth)}\\
r & \text{axisymmetric (\textbf{per radian}, } r \text{ the radial coordinate)}
\end{cases}
$$

This matches the code exactly: `CorrectParticleIntegrationWeightAxisymmetric`
(`src/AdjustParticleDiscretisation.FOR:183`) multiplies every material point's integration
weight by $r$ once at initialisation, and the internal force assembly scales the Gauss weight
by `Position` (`src/MPMDYNBTSig.FOR:191`). **All axisymmetric quantities in Anura3D are
per-radian** — multiply by $2\pi$ for a full solid of revolution. The rigid-body mass follows
the same convention (`RigidBodyMass`, `src/StructuralElements.for`).

Axisymmetric analysis adds the hoop strain and hoop stress:

$$
\varepsilon_{\theta\theta} = \frac{u_r}{r},
\qquad
(\nabla\!\cdot\!\boldsymbol\sigma)_r
= \frac{\partial \sigma_{rr}}{\partial r}
+ \frac{\partial \sigma_{rz}}{\partial z}
+ \frac{\sigma_{rr}-\sigma_{\theta\theta}}{r}
$$

which is why the internal force carries the extra term $\sigma_{\theta\theta} N_a / r$
(`src/MPMDYNBTSig.FOR:244`).

> **Note on the linear triangle.** For TRI3, $\mathbf{B}_a$ is *constant over the element*, so
> $\varepsilon_{rr}$, $\varepsilon_{zz}$, $\gamma_{rz}$ are element-constant. In axisymmetry
> $\varepsilon_{\theta\theta}=u_r/r$ is the **only** strain component that varies within an
> element. This single fact drives two things later: why the volumetric strain smoothing is
> mandatory (§8.2), and why the ghost-penalty face integral collapses to an algebraic
> expression (§4.3).

---

## 2. Governing equations

### 2.1 Strong form

In each sub-domain $\Omega^\pm$:

$$
\rho \frac{D\mathbf{v}}{Dt} = \nabla\!\cdot\!\boldsymbol\sigma + \rho\,\mathbf{g}
\qquad\text{in } \Omega^\pm
$$

with the usual boundary conditions on $\partial\Omega$, and on the interface $\Gamma$ the
**unilateral contact conditions**. Writing the penetration and the relative velocity as

$$
d = -\phi \quad (\text{contact active when } d>0),
\qquad
\dot{\mathbf{g}} = \mathbf{v} - \mathbf{v}_R,
\qquad
\dot g_n = \dot{\mathbf{g}}\cdot\mathbf{n},
\qquad
\dot{\mathbf{g}}_t = \dot{\mathbf{g}} - \dot g_n \mathbf{n}
$$

and the traction on the soil as $\mathbf{t} = \boldsymbol\sigma\mathbf{n} = t_n \mathbf{n} + \mathbf{t}_t$:

**Signorini (non-penetration):**

$$
d \le 0, \qquad t_n \le 0, \qquad t_n\, d = 0
$$

**Coulomb friction:**

$$
\|\mathbf{t}_t\| \le \mu |t_n|,
\qquad
\begin{cases}
\dot{\mathbf{g}}_t = \mathbf{0} & \text{if } \|\mathbf{t}_t\| < \mu|t_n| \quad(\text{stick})\\[2pt]
\mathbf{t}_t = -\mu|t_n| \dfrac{\dot{\mathbf{g}}_t}{\|\dot{\mathbf{g}}_t\|} & \text{otherwise} \quad(\text{slip})
\end{cases}
$$

### 2.2 Weak form with a jump across $\Gamma$

Find $\mathbf{v}$ such that for all admissible $\mathbf{w}$, **allowing $\mathbf{w}$ and
$\mathbf{v}$ to be discontinuous across $\Gamma$**:

$$
\underbrace{\sum_{i\in\{+,-\}}\int_{\Omega^i}\!\! \rho\,\mathbf{w}\cdot\mathbf{a}\;\varrho\,d\Omega}_{\text{inertia}}
+ \underbrace{\sum_{i}\int_{\Omega^i}\!\! \nabla^{s}\mathbf{w} : \boldsymbol\sigma\;\varrho\,d\Omega}_{\text{internal}}
=
\underbrace{\sum_{i}\int_{\Omega^i}\!\! \rho\,\mathbf{w}\cdot\mathbf{g}\;\varrho\,d\Omega
+ \int_{\Gamma_t}\!\! \mathbf{w}\cdot\bar{\mathbf{t}}\;\varrho\,d\Gamma}_{\text{external}}
+ \underbrace{\int_{\Gamma} \mathbf{w}\cdot\mathbf{t}\;\varrho\,d\Gamma}_{\text{interface}}
$$

**This is the whole point of the method.** The sum over $\Omega^+$ and $\Omega^-$ replaces a
single integral over $\Omega$, the test function is allowed to jump at $\Gamma$, and the
interface integral is the only coupling between the two sides.

In axisymmetry the internal term expands to include the hoop contribution:

$$
\nabla^{s}\mathbf{w} : \boldsymbol\sigma
= \frac{\partial w_r}{\partial r}\sigma_{rr}
+ \frac{\partial w_z}{\partial z}\sigma_{zz}
+ \left(\frac{\partial w_r}{\partial z}+\frac{\partial w_z}{\partial r}\right)\sigma_{rz}
+ \frac{w_r}{r}\,\sigma_{\theta\theta}
$$

and with $\varrho = r$ the last term contributes $w_r \sigma_{\theta\theta}$ — finite at the
axis, which is why the $r$ in $\varrho$ and the $r$ in the denominator must cancel *before*
evaluation (see the guard at `src/MPMDYNBTSig.FOR:235-245`).

### 2.3 Nitsche's form of the interface traction

Nitsche imposes the contact conditions weakly through a consistent term instead of a Lagrange
multiplier. Define the **augmented normal stress**

$$
P_\gamma(\mathbf{v}) \;=\; \sigma_n(\mathbf{v}) \;-\; \gamma\, d ,
\qquad
\gamma = \frac{\gamma_0 E}{h}
$$

Then the interface contribution becomes

$$
\boxed{\;
\int_{\Gamma}\mathbf{w}\cdot\mathbf{t}\;\varrho\,d\Gamma
\;\longrightarrow\;
\frac{1}{\gamma}\int_{\Gamma}\big[P_\gamma(\mathbf{v})\big]_{-}
\Big(\theta\,\sigma_n(\mathbf{w}) - \gamma\, w_n\Big)\;\varrho\,d\Gamma
\;}
$$

with $[x]_- = \min(x,0)$ and $\theta \in \{1,0,-1\}$ selecting the symmetric, incomplete and
skew-symmetric variants. The condition is satisfied automatically: in contact
($d=0$, $\sigma_n<0$) the bracket is active; when separated ($d<0$, $\sigma_n=0$) it vanishes.

**Relation to Anura3D's current law — this is the key to the staged plan.** Take $\theta = 0$
and drop the resolved stress $\sigma_n(\mathbf{v})$ from the bracket:

$$
\frac{1}{\gamma}\int_\Gamma [-\gamma d]_-\,(-\gamma w_n)\,\varrho\,d\Gamma
= \gamma\int_\Gamma d\, w_n \,\varrho\,d\Gamma
\qquad (d>0)
$$

which is **exactly a penalty spring on the penetration**. Anura3D's `LevelSetContact` computes

```fortran
Pmag = K_eff * k_of_d * (-1.0 * DistanceField(I)) - DampCoeff * vn
```

i.e. $\gamma \leftrightarrow K_{\text{eff}}\,k(d)$ with a log-barrier modulation $k(d)$, plus a
dashpot. **So the existing contact law is the $\theta=0$, stress-free limit of Nitsche.**
Phase 2 therefore *adds* the $\sigma_n$ consistency term rather than replacing the law, and
$\theta$ and the stress term can be switched on independently for comparison.

Friction enters the same way, with the tangential trial traction projected onto the Coulomb
disc of radius $\mu|[P_\gamma]_-|$ — structurally the same return mapping already implemented
in `LevelSetContact` (`Tstore` / `Ttrial` / `Tlimit`), so the existing elastic-predictor /
plastic-corrector code carries over.

### 2.4 Ghost-penalty stabilisation

When $\Gamma$ clips a sliver off an element, the physical volume supporting some nodal shape
functions becomes arbitrarily small. The nodal mass and stiffness contributions vanish, the
bilinear form loses coercivity, and the explicit update divides by a near-zero mass. The cure
is to penalise the jump of the normal gradient across the faces near the interface, which ties
the badly-supported dofs to well-supported neighbours:

$$
\boxed{\;
j(\mathbf{w},\mathbf{v}) \;=\;
\gamma_A\, h \sum_{F\in\mathfrak{F}_G} \int_F
\big[\![\, \nabla\mathbf{w}\cdot\mathbf{n}_F \,]\!\big] \cdot
\big[\![\, \nabla\mathbf{v}\cdot\mathbf{n}_F \,]\!\big] \; \varrho\, ds
\;}
$$

$\gamma_A$ carries **units of stress**; take $\gamma_A = \gamma_0^{\!A} E$ and calibrate
$\gamma_0^{\!A}$ (§8.1). The stabilised face set $\mathfrak{F}_G$ is built in two steps:

1. **Cut elements** — elements the interface passes through (§5.1).
2. **Stabilised faces** — faces shared between a cut element and a non-cut active element, or
   between two adjacent cut elements.

Following the MPM ghost-stabilisation literature, $j$ is integrated with **Gauss quadrature on
the grid, not with material points** — deliberately, because MP quadrature is exactly what is
unreliable inside a sliver.

---

## 3. The level set

### 3.1 Definition and orientation

For a bar chain $\Gamma = \bigcup_e \Gamma_e$ with segments $\Gamma_e = [\mathbf{X}_1^e,\mathbf{X}_2^e]$,
the signed distance at a point $\mathbf{x}$ is

$$
\phi(\mathbf{x}) = s(\mathbf{x})\,\min_e \operatorname{dist}(\mathbf{x},\Gamma_e),
\qquad
s(\mathbf{x}) = \operatorname{sign}\!\big((\mathbf{x}-\mathbf{x}^\star)\cdot \mathbf{n}_{e^\star}\big)
$$

where $e^\star$ is the closest segment, $\mathbf{x}^\star$ its closest point, and
$\mathbf{n}_{e^\star}$ its **unit** normal with a **consistent orientation along the chain**.

Distance to one segment, with $\mathbf{t}_e$ the unit tangent and $L_e$ the length:

$$
\tau = \frac{(\mathbf{x}-\mathbf{X}_1)\cdot\mathbf{t}_e}{L_e},
\qquad
\mathbf{x}^\star =
\begin{cases}
\mathbf{X}_1 & \tau \le 0\\
\mathbf{X}_1 + \tau L_e \mathbf{t}_e & 0<\tau<1 \quad(\text{face region})\\
\mathbf{X}_2 & \tau \ge 1 \quad(\text{end region})
\end{cases}
$$

### 3.2 Three corrections required in the existing code

The current per-MP field is **not** a signed distance and cannot be used as-is:

1. **`NormalsRigidBody` is scaled by `Thickness`, not unit**
   (`src/StructuralElements.for:347-354`). Consequently
   `D1 = abs(DotProduct(Xp-X1, Normal, 2))` (`:632`) evaluates to
   $\text{Thickness}\times\text{distance}$, yet in the same min-reduction it is compared
   against true Euclidean node distances from `Length(...)`. For $\text{Thickness}\ne 1$ the
   face-vs-end branch selection is biased. **The normal must be normalised.**
2. **The field is unsigned.** `DistanceField` is a minimum over unsigned distances, minus the
   particle radius. The sign is re-derived ad hoc by dot products in three separate places
   (`:682-694`, `:1033-1046`, `:1251-1256`). A single consistent orientation per body is
   required, computed once.
3. **There is no nodal level set at all** — only a per-MP one. The compatibility test of §4.1
   is a *node versus particle* test, so $\phi_a$ must be added.

### 3.3 Precondition: the chain must not end inside the mesh

A signed distance is single-valued only for a chain that is **closed**, or that **terminates on
the domain boundary or on the symmetry axis**. A chain ending strictly inside the mesh creates
a branch cut emanating from the free end, across which the sign flips spuriously and elements
are wrongly reported as cut.

For the axisymmetric cone this is satisfied — the chain runs from the axis ($r=0$) to the outer
or upper boundary. **This should be asserted at input** rather than assumed. Handling genuinely
interior free ends (a crack tip) requires branch enrichment and is out of scope.

---

## 4. Material point discretisation

### 4.1 The compatibility condition

This is the CPIC construction. Each node carries a side tag, each material point carries one,
and a particle–node pair is **compatible** when they are on the same side:

$$
\chi_{pa} =
\begin{cases}
1 & \operatorname{sign}\phi_p = \operatorname{sign}\phi_a \quad\text{(or either is outside the band)}\\
0 & \text{otherwise}
\end{cases}
$$

For several rigid bodies the test must hold **for every body whose band contains both** the
particle and the node — this is the "colored" part of the colored distance field, stored as a
small bitmask per node and per particle.

> **Why no extra degrees of freedom are needed.** Each *node* has a definite sign of $\phi$, so
> it belongs to exactly one side. The discontinuity lives strictly *inside* cut elements, never
> at a node. The two-sidedness is produced by **masking the transfers**, not by storing two
> values per node. `TotalVelocitySoil` stays a single array, `ReducedDof` is unchanged, the
> solver is unchanged, and the existing multi-entity contact machinery is untouched.

### 4.2 Masked transfers

**Particle to grid — a scatter, so no partition of unity is required.** An incompatible node
simply receives nothing from that particle, which is precisely the intent: a node's mass is the
mass of its own side only.

$$
m_a = \sum_p m_p\, N_a(\mathbf{x}_p)\, \chi_{pa},
\qquad
(m\mathbf{v})_a = \sum_p m_p\, \mathbf{v}_p\, N_a(\mathbf{x}_p)\, \chi_{pa}
$$

$$
\mathbf{f}^{int}_a = \sum_p \mathbf{B}_a^{\!\top}\!(\mathbf{x}_p)\,\boldsymbol\sigma_p\; w_p\, \chi_{pa},
\qquad w_p = \text{`Particles(p)\%IntegrationWeight'} \;(\propto r \text{ in axisymmetry})
$$

**Grid to particle — a gather, so $\sum_a N_a = 1$ must be preserved.** Do *not* drop the
incompatible terms; substitute the rigid-body velocity at that node:

$$
\boxed{\;
\mathbf{v}_p = \sum_a N_a(\mathbf{x}_p)\Big[\chi_{pa}\,\mathbf{v}_a + (1-\chi_{pa})\,\mathbf{v}_R(\mathbf{x}_a)\Big]
\;}
$$

Partition of unity holds exactly, and the physics is right: a particle hugging the bar is
bounded by the bar, so the nodes on the far side contribute the bar's own motion rather than
soil motion that never reached it.

**Degenerate case — the small cut.** If a bar clips a corner so that all three nodes lie on one
side, then $\sum_a \chi_{pa}N_a = 0$ for the minority particles and the formula above gives
$\mathbf{v}_p = \mathbf{v}_R$: those particles move rigidly with the bar. Their contribution to
nodal mass also vanishes. This is exactly the ill-conditioning that ghost penalty exists to
cure — **the two phases meet here**, and a nodal mass floor (§8.3) is mandatory before Phase 1
can be trusted.

### 4.3 Ghost penalty for element-constant gradients

Because TRI3 gradients are element-constant, $\nabla\mathbf{v}|_K = \sum_b \mathbf{v}_b \otimes \mathbf{B}_b|_K$
is constant on each element, so the jump is constant on the face and the integral is exact
without quadrature. For a face $F$ between elements $K_1$ and $K_2$, define the scalar

$$
\beta_b^F = \big(\mathbf{B}_b|_{K_1} - \mathbf{B}_b|_{K_2}\big)\cdot \mathbf{n}_F
\qquad
(\mathbf{B}_b|_K = \mathbf{0} \text{ if node } b \notin K)
$$

so that the gradient jump is $[\![\nabla\mathbf{v}\cdot\mathbf{n}_F]\!] = \sum_b \beta_b^F\,\mathbf{v}_b$.
The stabilisation and the resulting nodal force are then purely algebraic:

$$
j = \gamma_A\, h \sum_{F}\; \varrho_F |F| \Big(\sum_a \beta_a^F \mathbf{w}_a\Big)\cdot\Big(\sum_b \beta_b^F \mathbf{v}_b\Big),
\qquad
\boxed{\;\mathbf{f}^{ghost}_a = \gamma_A\, h \sum_{F} \varrho_F |F|\, \beta_a^F \sum_b \beta_b^F \mathbf{v}_b \;}
$$

with $|F|$ the edge length and $\varrho_F = 1$ (plane strain) or $\bar r_F$, the edge midpoint
radius (axisymmetric, per radian). **No face quadrature loop is required** — this is the single
biggest implementation saving available.

$\mathbf{f}^{ghost}$ is subtracted alongside the internal force at
`src/LagrangianPhase.FOR:485`.

### 4.4 Critical time step

The ghost term adds stiffness, so it *reduces* the stable explicit step. Since
$\beta \sim 1/h$ and $|F| \sim h$, the added stiffness scales as
$\gamma_A h \cdot h \cdot h^{-2} = \gamma_A$, to be compared against the element stiffness
$\sim E$. So

$$
\Delta t_{crit} \;\lesssim\; \mathrm{Cr}\,\frac{h}{c},
\qquad
c = \sqrt{\frac{K + \tfrac{4}{3}G + c_\ast \gamma_A}{\rho}}
$$

with $c_\ast$ an $O(1)$ constant to be measured, not assumed. With
$\gamma_A = \gamma_0^{\!A} E$ this is a penalty of roughly
$\sqrt{1+\gamma_0^{\!A}}$ on the time step. **This must be fed into
`CalculateCriticalTimeStep` (`src/MPMDynamicExplicit.FOR:154`)** — otherwise enabling ghost
penalty silently destabilises the run it was meant to stabilise.

---

## 5. Cut detection and the narrow band

### 5.1 The cut test

The bars **are** the cut geometry — no implicit surface reconstruction, no marching squares,
no segment/triangle intersection code. Evaluate the level set at the **nodes, preserving sign**;
then

$$
\boxed{\;K \text{ is \textbf{cut}} \iff \min_{a\in K}\phi_a \;<\; 0 \;<\; \max_{a\in K}\phi_a\;}
$$

The same nodal signs serve as the CPIC side tag of §4.1, so one field supplies both.

### 5.2 Why a narrow band is necessary

The current search is brute force — every particle against every bar,
$O(N_{MP}\!\times\!N_{bar})$ (`src/StructuralElements.for:611-669`) — with no spatial culling
anywhere in the code. Worse, **the same loop is duplicated verbatim three times**
(`DistanceToRigidBody:611-669`, `LevelSetContact:972-1017`,
`DisplacementCorrection:1344-1390`) and `LevelSetContact` runs it **every time step**
(`src/MPMDynamicExplicit.FOR:103`). Computing $\phi$ at every node instead would merely move
the same cost.

### 5.3 Construction by breadth-first sweep

1. **Seed.** For each bar node (from `ThinRigidCoordinates`), find its containing element with
   `DetermineGlobPosElement` (`src/ElemConnections.FOR:2428`), seeded with **last step's
   element** for that bar node. With temporal coherence this is 0–1 hops.
2. **Sweep.** Breadth-first over `ElementAdjacencies` / `GetAdjacentElement`
   (`src/ElemConnections.FOR:535`). For each newly reached element compute $\phi_a$ at its
   nodes against a **propagated candidate bar list** — the bars carried by the element it was
   reached from. Stop expanding through an element once $\min_a|\phi_a| > R_{band}$.
3. **Classify.** `ElementBandDepth`: *cut* (sign change), *band* ($|\phi|\le R_{band}$), *far*.

Because $\phi$ is a true distance to the bar set, the stopping criterion makes the sweep follow
the chain automatically; **one seed per connected chain suffices**.

### 5.4 Band radius — derived, not guessed

The band must be wide enough that no contact can develop between two rebuilds. Between
rebuilds a particle can approach the bar by its own radius, and the bar itself moves:

$$
\boxed{\;
R_{band} \;\ge\; \max_p \big(\text{ParticleRadius}_p\big)
\;+\; \|\mathbf{v}_{bar}\|\,\Delta t\, n_{rebuild}
\;+\; \delta_{safety}
\;}
$$

With rebuild every step ($n_{rebuild}=1$) the middle term is one particle displacement and is
negligible; the criterion is then dominated by the particle radius. Express $R_{band}$ in units
of `ElementLMin` so it is mesh-independent, and **assert the inequality at runtime** — a band
that is silently too narrow produces missed contacts, not a crash.

This is the requirement that the influence region be "larger than a single element width":
$R_{band} \gtrsim$ `ElementLMin` follows from the particle radius being an appreciable fraction
of the element size.

### 5.5 Data structures

New module `ModCutCellTopology` (new file, added to `src/VS/Anura3D.vfproj`):

```fortran
real(REAL_TYPE),    allocatable :: NodalPhi(:)          ! signed level set, valid in the band
integer(INTEGER_TYPE), allocatable :: NodalSide(:)      ! +1 / -1 / 0 (far)  -- the CPIC tag
integer(INTEGER_TYPE), allocatable :: ElementBandDepth(:)   ! cut / band / far
integer(INTEGER_TYPE), allocatable :: BandElementList(:)
integer(INTEGER_TYPE)              :: NBandElements
integer(INTEGER_TYPE), allocatable :: ElementBarOffsets(:), ElementBarList(:)  ! CSR
integer(INTEGER_TYPE), allocatable :: BarSeedElement(:)  ! temporal coherence for the seed
```

Rebuilt once per step, between `RigidBodyLagrangianPhase` and `LevelSetContact`
(`src/MPMDynamicExplicit.FOR:102-103`).

### 5.6 Rewiring the per-MP field

Loop `BandElementList`, take its particles via `NPartEle` / `GetParticleIndex`, and test each
only against `ElementBarList` for its own element. Particles outside the band get
$\phi_p = +\infty$ and a null affinity. **Expected speedup 2–3 orders of magnitude** for a CPT
model ($\sim\!10^5$ particles, $\sim\!10^2$ bars): cost becomes proportional to the particles
actually near the cone rather than to the whole domain. Collapse the three duplicate searches
into this one routine in the same pass.

Housekeeping note: `AffinityArray` is declared `real(INTEGER_TYPE)` — a single-precision
**real** used directly as an array subscript (`src/MPMData.FOR:104`). Extending it is the
moment to make it `integer(INTEGER_TYPE)`.

---

## 6. Pseudo-algorithm for one time step

Phase 1 lines are unmarked; **Phase 2 additions are marked `[P2]`** and can be omitted entirely
without affecting the rest.

```
--- RIGID BODY ---------------------------------------------------------------
 1  RigidBodyLagrangianPhase()                       [StructuralElements.for]
       Mass  = RigidBodyMass(J)                       (Mg, per radian if axisym)
       a_R   = (Fsum + Weight + Fprescribed) / Mass
       v_R  += a_R dt ;  move ThinRigidCoordinates ;  refresh unit normals

--- CUT TOPOLOGY  (NEW) ------------------------------------------------------
 2  BuildCutTopology()                               [ModCutCellTopology]
       for each bar node: seed element <- DetermineGlobPosElement(seed = last step)
       BFS over ElementAdjacencies:
           phi_a  <- signed distance to candidate bars      (Sec. 3.1)
           stop expanding when min_a |phi_a| > R_band       (Sec. 5.4)
       NodalSide(a) <- sign(phi_a)
       ElementBandDepth(K) <- cut  if  min phi_a < 0 < max phi_a
       assert R_band >= max ParticleRadius + |v_bar| dt

 3  UpdateLevelSet()                  (replaces 3 duplicated searches)
       for K in BandElementList, for p in K:
           phi_p <- signed distance to ElementBarList(K)
           chi(p,a) <- [ sign(phi_p) == NodalSide(a) ]      (Sec. 4.1)

--- INTERFACE TRACTION -------------------------------------------------------
 4  InterfaceTraction()                              [was LevelSetContact]
       Phase 1: keep the existing penalty law, evaluated per MP
                Pmag = K_eff k(d) d - c_damp vn ;  Coulomb return mapping
  [P2]          integrate over the CUT SEGMENT instead of per MP
  [P2]          add the Nitsche consistency term  (1/g)[P_g(v)]_- theta sigma_n(w)   (Sec. 2.3)
       accumulate Fsum, Tsum on the rigid body      (two-way coupling)
       ContactForceArray(p,:) = -Force

--- LAGRANGIAN PHASE ---------------------------------------------------------
 5  MapMomentumAndMassP2N()                          [LagrangianPhase.FOR:1244]
       m_a      += m_p N_a(x_p) * chi(p,a)                  <-- MASKED
       (mv)_a   += m_p v_p N_a(x_p) * chi(p,a)              <-- MASKED

 6  v_a = (mv)_a / m_a          with a MASS FLOOR, not /= 0   (Sec. 8.3)

 7  MPMDYNBTSig()                                    [MPMDYNBTSig.FOR:61]
       f_int_a += B_a^T sigma_p w_p * chi(p,a)              <-- MASKED
       axisym: + sigma_tt N_a / r      (guarded at r = 0)

  [P2] GhostPenaltyForce()                            (Sec. 4.3)
  [P2]     for F in stabilised faces:
  [P2]         beta_b = (B_b|K1 - B_b|K2) . n_F
  [P2]         f_ghost_a += gamma_A h |F| rho_F beta_a sum_b beta_b v_b

 8  GetNodalContactForce(Fcontact)                   [StructuralElements.for:1427]
       Fcontact_a += ContactForceArray(p,:) N_a(x_p) * chi(p,a)   <-- MASKED

 9  RateofMomentum = (ExtLoad + GravityLoad - IntLoad + Fcontact
  [P2]                - f_ghost ) / ScalingMassFactor  [LagrangianPhase.FOR:485]
10  a_a = RateofMomentum / m_a        with the same MASS FLOOR
11  v_a += a_a dt ;  apply PBoundary fixities

--- CONVECTIVE PHASE ---------------------------------------------------------
12  UpdateParticleVelocityAndMapMomentum()           [MPMDYNConvPhase.FOR:944]
       v_p = sum_a N_a [ chi(p,a) v_a + (1-chi(p,a)) v_R(x_a) ]   <-- MASKED + SUBSTITUTED
       (APIC/FLIP blend applies to the masked gather, not around it)
       if sum_a chi(p,a) N_a == 0  ->  v_p = v_R      (degenerate small cut)

13  strain update:  Get_Strain per MP, or SmoothenStrains
       CUT-AWARE: the nodal patch average must not mix sides   (Sec. 8.2)

14  MPMDYNGetSig()  ->  constitutive update
15  UpdateParticlePos, UpdateParticleHouseKeeping, SetActiveElement
16  CalculateCriticalTimeStep()
  [P2]   include the ghost stiffness in the wave speed          (Sec. 4.4)
```

---

## 7. Difficulty assessment and staging

| Item | Difficulty | Estimate |
| --- | --- | --- |
| Nodal level set + cut/band detection (§3, §5) | Easy–Moderate | 1–2 weeks |
| Prerequisite fixes (§8.3, §8.4) | Easy | a few days |
| **Phase 1** — CPIC compatibility (§4.1, §4.2) | Moderate | 4–6 weeks |
| Cut-aware volumetric strain smoothing (§8.2) | Moderate–Hard | 2–3 weeks |
| **Phase 2** — ghost penalty (§4.3, §4.4) | Hard | 4–6 weeks |
| **Phase 2** — Nitsche consistency term (§2.3) | Hard | 4–8 weeks |

**Why Phase 1 is only moderate.** Anura3D's transfer loops already expose the local node index,
so the mask drops straight in at six places:

| Site | File | Operation |
| --- | --- | --- |
| P2G mass | `src/LagrangianPhase.FOR:1301` | `LumMassX(:) += ParticleMass * PartShape(:)` |
| P2G momentum | `src/LagrangianPhase.FOR:1305` | `Mom(I,:) += ParticleMass * PartShape(:) * v(I)` |
| G2P acceleration | `src/MPMDYNConvPhase.FOR:1039` | `+= ParticleShape(INode) * NodAcc(...)` |
| G2P APIC velocity | `src/MPMDYNConvPhase.FOR:1063` | `+= ParticleShape(INode) * NodalVelocity(...)` |
| G2P momentum remap | `src/MPMDYNConvPhase.FOR:1085` | `Momentum(...) += Mass * ParticleShape * v(I)` |
| Contact force scatter | `src/StructuralElements.for:1459-1460` | `Fcontact += ContactForceArray * PartShape(INode)` |

No new degrees of freedom, no change to `ReducedDof` (`src/MeshInfo.FOR:1994-2025`), no change
to `Counters%NEntity`, no change to the solver, output or restart. The Bardenhagen multi-entity
contact algorithm is left untouched and remains the mutually-exclusive alternative it is today.

---

## 8. Open questions and prerequisites

### 8.1 Calibration of $\gamma_A$ and $\gamma_0$

$\gamma_A$ has units of stress and $\gamma = \gamma_0 E/h$ has units of stress per length.
Both are accuracy *and* stability parameters here (§4.4). Too small and the small-cut
conditioning is not fixed; too large and the time step collapses. Needs a parametric study on
the sliver benchmark of §9.

### 8.2 Do ghost penalty and volumetric strain smoothing compose, or fight?

Both are stabilisations acting on the same elements, and the axisymmetric cone model
**depends on the volumetric strain smoothing to run at all** — TRI3 in axisymmetry is severely
volumetrically locking because $\varepsilon_{\theta\theta}$ is the only non-constant strain
component, so the incompressibility constraint is enforced closer to pointwise than in plane
strain. Two specific concerns:

- The nodal patch average in `GetEnhancedVolumetricStrain` (`src/MPMStrainSmoothing.FOR:81`)
  accumulates `StrainSmooth(node,1) / StrainSmooth(node,2)` over **all** elements touching a
  node. A node on the $+$ side must accumulate only from $+$-side contributions, otherwise the
  smoothing reconnects across the cut everything the masking just disconnected. **This is on
  the critical path for Phase 1, not an afterthought.**
- Ghost penalty penalises the gradient jump across faces; volumetric smoothing averages the
  volumetric strain across the same patch. Whether the combination is over-stiff has to be
  measured.

### 8.3 A nodal mass floor is mandatory

The only guard today is an exact-zero test, `if (LumpedMassDry(IDOF,J)/=0)`, at
`src/LagrangianPhase.FOR:1423, 1436, 1469, 1482`. Once masking removes contributions, nodal
masses become small but non-zero and the division amplifies noise. A relative floor
(e.g. a fraction of the mean nodal mass of the patch) must be added at all four sites **before**
Phase 1 is trusted.

### 8.4 Two pre-existing defects to clear first

- **`src/MPMDYNConvPhase.FOR:1015`** —
  `NodalVelocity(I,:) = TotalVelocitySoil(IDof(I,:), 1)` hard-codes entity 1 while the FLIP
  branch at `:1039` uses `iEntity`. Harmless while `NEntity == 1`, but this is the line that
  must carry the mask, so it gets rewritten anyway.
- **`src/StructuralElements.for:1441`** — `GetNodalContactForce` declares
  `Fcontact(Counters%N, 2)` while the caller declares
  `FContact(Counters%N, Counters%nEntity)` (`src/LagrangianPhase.FOR:84`). With contact off
  (`NEntity == 1` — **the current models**) the callee writes past the end of the array.
  **This is a live bug in the runs being done today** and should be fixed independently of this
  project.

### 8.5 Friction in Nitsche form

The Coulomb return mapping already implemented in `LevelSetContact` is an
elastic-predictor/plastic-corrector on a tangential penalty spring. Its Nitsche counterpart
projects the augmented tangential traction onto the Coulomb disc. The two coincide when the
stress consistency term is dropped, so the migration is incremental — but the friction
parameter may need recalibration once $\sigma_n$ enters.

---

## 9. Verification plan

| Test | Purpose | Pass criterion |
| --- | --- | --- |
| **Zero-cut regression** | the band machinery must be inert when unused | with no rigid body, results bit-identical to today |
| **Level-set unit test** | §3, §5.1 | one straight bar over a known TRI3 patch: cut set matches the analytical one; band symmetric about the bar |
| **Band radius assertion** | §5.4 | runtime assert fires when $R_{band}$ is set below the derived bound |
| **Compatibility unit test** | §4.1, §4.2 | for an element straddling the bar, the mask selects the expected nodes and the G2P weights still sum to 1 |
| **Cost measurement** | §5.6 | wall-clock of `LevelSetContact` before/after; `DistanceField` unchanged for every particle inside the band |
| **Rigid-body statics** | no regression in the mass/weight fix | startup log line `Rigid body 1 (...) axisymmetric: mass = ... Mg/rad` unchanged |
| **Phase 1 physical target** | the reason for the project | axisymmetric cone: material points behind the shoulder separate from the shaft while those ahead of the tip move with it; compare shaft friction and tip resistance against the current single-field result |
| **Phase 2 sliver test** | §2.4 | place a bar to clip an arbitrarily thin sliver: unstable without ghost penalty, stable with it |
| **Phase 2 patch test** | Nitsche consistency | uniform stress state transmitted exactly across a cut element for each $\theta \in \{1,0,-1\}$ |

---

## 10. References

**CutFEM and ghost penalty**
- Burman & Hansbo et al., [*CutFEM based on extended finite element spaces*](https://link.springer.com/article/10.1007/s00211-022-01313-z) ([preprint](https://arxiv.org/pdf/2101.10052))
- [*On the design of locking-free ghost penalty stabilization*](https://link.springer.com/article/10.1007/s00211-025-01502-6)
- [*Projected gradient stabilization of sharp and diffuse interface unfitted Nitsche methods*](https://arxiv.org/pdf/2501.16594)

**Nitsche for contact**
- Chouly & Hild, [*Nitsche-based FEM for contact with Coulomb friction*](https://link.springer.com/chapter/10.1007/978-3-319-96415-7_79)
- Chouly, Mlika & Renard, [*An unbiased Nitsche's approximation of the frictional contact between two elastic structures*](https://link.springer.com/article/10.1007/s00211-018-0950-x)
- [*Nitsche method for contact with Coulomb friction: existence results, static and dynamic*](https://www.sciencedirect.com/science/article/abs/pii/S0377042722002746)
- Mlika, Renard & Chouly, [*An unbiased Nitsche's formulation of large deformation frictional contact and self-contact*](https://www.researchgate.net/publication/318593789_An_unbiased_Nitsche's_formulation_of_large_deformation_frictional_contact_and_self-contact)
- Mlika, [*Nitsche method for frictional contact and self-contact* (thesis)](https://theses.hal.science/tel-02067118)

**Ghost stabilisation inside MPM**
- Coombs, [*Ghost stabilisation of the Material Point Method for stable quasi-static and dynamic analysis of large deformation problems*](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.7332), IJNME 2023 ([preprint](https://arxiv.org/pdf/2209.10955))
- Pretti, Coombs et al., [*A stable poro-mechanical formulation for MPM leveraging overlapping meshes and multi-field ghost penalisation*](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.7630) ([preprint](https://arxiv.org/html/2405.12814v1))

**Velocity discontinuities in MPM**
- Hu et al., [*A moving least squares material point method with displacement discontinuity and two-way rigid body coupling*](https://dl.acm.org/doi/10.1145/3197517.3201293), SIGGRAPH 2018 — the CPIC algorithm ([project page](https://yuanming.taichi.graphics/publication/2018-mlsmpm/), [reference implementation](https://github.com/Zhijie-YU/MLS_MPM_CPIC))
- Homel & Herbold, [*Modeling strong discontinuities in the material point method using a single velocity field*](https://www.sciencedirect.com/science/article/abs/pii/S0045782518305607), CMAME 2019

**Application context**
- [*Investigation of the Material Point Method in the simulation of Cone Penetration Tests in dry sand*](https://www.sciencedirect.com/science/article/abs/pii/S0266352X20304869)
- [*Two-phase Material Point Method applied to the study of cone penetration*](https://www.sciencedirect.com/science/article/abs/pii/S0266352X16300477)
- [*Relating cone penetration resistance to sand state through MPM large-deformation modelling*](https://www.researchgate.net/publication/361008245_Relating_cone_penetration_resistance_to_sand_state_through_material_point_method_large-deformation_modelling)
