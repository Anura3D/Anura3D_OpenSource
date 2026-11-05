# Gasification of stress-free material points in Anura3D

**Status:** implemented (`src/Gasification.FOR`), explicit dynamic, single-phase. See §9 for what
was built and where it differs from the design below.
**Scope:** 2D plane strain and 2D axisymmetric on TRI3 meshes; 3D TET4 works the same way.
Explicit dynamic MPM-MIXED, dry single-phase material (`NFORMULATION = 1`, one phase), with or
without Bardenhagen contact and the rigid body. Saturated and unsaturated formulations are out of
scope for the first version.
**Target application:** ejecta and loosened soil thrown out by a cone or wall impact.

---

## 1. Concept

An impact throws soil out of the main body. These material points (MPs) quickly lose all
confinement: their stress goes to zero and the density drops. As a continuum they are a problem:

* each sparse MP activates an element, so many nearly empty elements carry stiffness;
* those elements are cut elements with tiny nodal masses, the source of the instability treated by
  ghost stabilisation (`Ghost_Stabilisation_Anura3D.md`);
* the constitutive model is driven at near-zero stress, where many soil models are badly conditioned.

**Gasification** switches such MPs to a **gas** state, in which they behave as discrete particles.
When they land and recompact, they switch back to **continuum** and the constitutive model is
reactivated.

Each MP is in one of two states:

| | Continuum (normal) | Gas |
| --- | --- | --- |
| Maps mass and momentum to the grid | yes | **yes** (unchanged) |
| Stress $\boldsymbol\sigma_p$ | from the constitutive model | **0** |
| Internal force contribution | $V_p\boldsymbol\sigma_p\nabla N_a$ | **none** |
| Constitutive model called | yes | **no** |
| Volume $V_p$ tracked | yes | **yes**, needed for reactivation |
| Gravity, contact with the body | yes | yes |
| Wave-speed time-step limit | yes | **no**; a kinematic limit applies instead |

Because a gas MP stays on the grid with its mass and momentum (user decision), total momentum is
conserved exactly. The Bardenhagen contact still sees gas MPs as soil and keeps them out of the cone.
In nodes shared with the continuum, gas MPs act as passive inertia that the continuum carries along.

---

## 2. Notation and measures

| Symbol | Meaning | Anura3D variable |
| --- | --- | --- |
| $V_p$, $m_p$ | MP volume and mass | `Particles(p)%IntegrationWeight`, `MassArray(p)` |
| $V_{p,0}$ | reference volume (initial, or volume at last reactivation) | new `Particles(p)%V0Ref` |
| $\chi_p=V_{p,0}/V_p$ | density ratio $\rho_p/\rho_{p,0}$ | computed |
| $n_p$, $n_{p,0}$ | porosity, reference porosity | `Particles(p)%Porosity`, `%InitialPorosity` |
| $p=\tfrac13\operatorname{tr}\boldsymbol\sigma'$ | mean effective stress, **tension positive** (Anura convention, consistent with `DynUpdateParticlePorosity`) | `SigmaEffArray` |
| $q=\sqrt{\tfrac32\,\mathbf s:\mathbf s}$ | deviatoric stress | |
| $\dot\varepsilon_v$ | volumetric strain rate, $\operatorname{tr}\Delta\boldsymbol\varepsilon/\Delta t$, positive in dilation | `GetEpsStepI` |
| $\phi_K$ | continuum fill of element $K$: $\sum_{p\in K,\,\text{continuum}}V_p/\lvert K\rvert$ | `ElementSpace` |
| $p_{\text{ref}}$ | reference pressure (user input, e.g. 1 kPa) | new input |

**Porosity from density.** The solid volume is conserved, $(1-n_{p,0})V_{p,0}=(1-n_p)V_p$, so

$$
n_p=1-(1-n_{p,0})\,\chi_p,\qquad e_p=\frac{n_p}{1-n_p}.
$$

This gives the porosity and void ratio for density-dependent models at reactivation, even when
`APPLY_POROSITY_UPDATE` is off.

---

## 3. Continuum → gas

An MP becomes gas when **all** of the following hold for $N_{\text{hold}}$ consecutive time steps.
The counter is reset as soon as any condition fails.

| # | Condition | Meaning |
| --- | --- | --- |
| G1 | $\lvert p\rvert\le\varepsilon_\sigma p_{\text{ref}}$ **or** $p>0$ | no confinement (soil carries no tension) |
| G2 | $q\le\varepsilon_\sigma p_{\text{ref}}$ | no shear stress |
| G3 | $\chi_p\le\chi_{\text{off}}$ **or** $\phi_K\le\phi_{\text{off}}$ | loosened, or sitting in a sparse region |
| G4 | $\dot\varepsilon_v\ge0$ | not being compressed |
| G5 | the material allows gasification (§7.2), the MP is not rigid body, and $K$ holds no rigid-body MPs or contact interface nodes | never gasify at the cone surface |

The hysteresis (G3 against R1 in §5) and the hold counter stop MPs from switching back and forth.

**On switching to gas:**

* store nothing new: mass, velocity, position and volume are kept, so momentum is unchanged;
* set $\boldsymbol\sigma_p=0$ and `IsBoundaryParticle = .true.`;
* freeze the constitutive state variables (they are reset on reactivation, §5.2).

---

## 4. Gas-state numerics

The gas state must not create new instabilities. Each item below names the routine it touches.

1. **Mass and momentum mapping: unchanged.**
   `MapMomentumAndMassP2N` (`LagrangianPhase.FOR:1244`) and `UpdateParticleVelocityAndMapMomentum`
   (`MPMDYNConvPhase.FOR:946`) treat gas MPs as normal soil MPs, so momentum is conserved exactly.
2. **Internal force: skip.**
   In `MPMDYNBTSig` (`MPMDYNBTSig.FOR`), gas MPs contribute no $V_p\boldsymbol\sigma_p\nabla N_a$ and
   no bulk viscosity load. In an element with only gas MPs, the internal force is zero and nodal
   accelerations are exactly $\mathbf g$ plus contact, so the flight is ballistic without grid noise.
3. **Integration type.** Any element that holds a gas MP is set to `IsParticleIntegration = .true.`
   in `CheckFillingOfElements` (`MPMConvPhase.FOR:952`). This keeps Gauss-point integration from
   treating gas mass and volume as stressed material. Gas MPs are not counted in the fill used there,
   so the continuum part of the element is classified on its own.
4. **Stress mapping: skip.**
   `AssignStressesToParticles` (`MPMConvPhase.FOR:1025`) must not copy Gauss-point stresses onto gas
   MPs, and `MPMDYNGetSig` must not call `StressSolid` (`MPMDYNStresses.FOR:274`, `:288`) for them.
5. **Volume.** The strain increment is still computed (`UpdateParticleStrains`) and $V_p$ is updated
   from it as usual, so $\chi_p$ follows expansion and recompaction. In very sparse regions the grid
   divergence is noisy. Reactivation therefore also uses the element fill (R2 in §5).
6. **Local and bulk damping.** Homogeneous local damping is proportional to the out-of-balance force.
   For gas nodes this force is gravity, and damping it would slow free flight. Gas mass is excluded in
   `MapDampingFromParticles` (`LagrangianPhase.FOR:2085`), so the damping coefficient of nodes with
   only gas mass is zero. Bulk viscosity is also skipped for gas MPs.
7. **Time step.** `CalculateCriticalTimeStep` (`MPMDYNConvPhase.FOR:2882`) skips the wave speed of
   gas MPs, as it already does for rigid-body MPs. It adds a kinematic limit so a gas MP cannot jump
   more than a fraction of an element per step:
   $$\Delta t\le C\,\frac{h_K}{\lvert\mathbf v_p\rvert+\epsilon}.$$
   Elements holding only gas MPs no longer constrain the step through $c$.
8. **Velocity transfer.** A gas MP uses the same FLIP update. As an option, the gas velocity is
   blended with the grid (PIC) velocity, $\mathbf v_p\leftarrow\alpha_g\mathbf v_p^{\text{FLIP}}+(1-\alpha_g)\mathbf v_p^{\text{PIC}}$,
   to damp noise in crowded gas regions. Here $\alpha_g=1$ is pure FLIP. The PIC velocity must be taken
   from the soil entity (`TotalVelocitySoil(:, SOFT_ENTITY)`).
9. **Ghost stabilisation.** Elements holding only gas MPs are excluded from ghost faces. Mixed
   elements are classified by their continuum fill.
10. **Gas mass in shared nodes.** An optional cap on the gas fraction of a node's mass is documented
    for later. It is not needed while gas MPs have no stiffness.

---

## 5. Gas → continuum (reactivation)

### 5.1 Conditions

A gas MP becomes continuum when either condition holds for $N_{\text{hold}}$ consecutive steps:

| # | Condition | Meaning |
| --- | --- | --- |
| R1 | $\chi_p\ge\chi_{\text{on}}$, with $\chi_{\text{on}}>\chi_{\text{off}}$ | the MP has recompacted |
| R2 | $\phi_K\ge\phi_{\text{on}}$ **and** $K$ holds continuum MPs **and** $\dot\varepsilon_v<0$ | it has landed in a filled region and is being compressed |

### 5.2 Re-initialisation

The aim is a stress-free "fresh deposit" that is consistent with the current density, with no change
in momentum.

1. **Kinematics:** keep $m_p$, $\mathbf v_p$, $\mathbf x_p$ and $V_p$. Momentum and kinetic energy are
   unchanged.
2. **Stress:** $\boldsymbol\sigma_p=0$. No stored elastic energy is introduced, so there is no shock.
3. **Density state:** $n_p=1-(1-n_{p,0})\chi_p$, written to `Particles(p)%Porosity`. Then reset the
   reference, $V_{p,0}\leftarrow V_p$ and $n_{p,0}\leftarrow n_p$, so later checks measure the change
   from the deposit.
4. **State variables:**
   * Mohr–Coulomb and linear elastic: reset accumulated plastic strain and any hardening variables
     to their initial values.
   * Density-dependent models: set the void ratio $e_p$ from step 3. For MCC, use the preconsolidation
     pressure on the normal compression line,
     $p_c=\exp\!\big((e_N-e_p)/\lambda\big)$, clipped at $p_{c,\min}$.
   * External (DLL) models: restore the MP's initial state-variable vector. A copy is stored at
     initialisation in a new array `ESMstatev0Array`, next to `ESMstatevArray` (`MPMData.FOR:140`).
     The void ratio is overwritten when the model exposes it.
5. **Stiffness ramp:** for the next $N_{\text{ramp}}$ steps, pass a scaled strain increment to the
   constitutive model,
   $$\Delta\boldsymbol\varepsilon_p\leftarrow\lambda_r\,\Delta\boldsymbol\varepsilon_p,\qquad
   \lambda_r=\min\!\big(1,\;k/N_{\text{ramp}}\big),$$
   where $k$ is the number of steps since reactivation. This works for every model, including
   external DLLs. It avoids a stiffness jump and keeps the time step from collapsing. Volume and
   porosity still use the full increment.
6. **Flags:** `IsGas = .false.`, `GasHoldCounter = 0`, and `IsBoundaryParticle = .true.` (the
   element is treated as partially filled until the next filling check).

---

## 6. Parameters

| Parameter | Meaning | Proposed start | Notes |
| --- | --- | --- | --- |
| $\varepsilon_\sigma$ | stress tolerance relative to $p_{\text{ref}}$ | 0.01 | |
| $p_{\text{ref}}$ | reference pressure | 1 kPa | of the order of the self-weight stress at one MP depth |
| $\chi_{\text{off}}$ | density ratio to gasify | 0.60 | about 40 % volume increase |
| $\chi_{\text{on}}$ | density ratio to reactivate | 0.80 | must be $>\chi_{\text{off}}$ |
| $\phi_{\text{off}}$ | element fill to gasify | 0.10 | |
| $\phi_{\text{on}}$ | element fill to reactivate | 0.40 | must be $>\phi_{\text{off}}$ |
| $N_{\text{hold}}$ | consecutive steps before switching | 5 | |
| $N_{\text{ramp}}$ | stiffness ramp length | 50 | |
| $\alpha_g$ | FLIP fraction for gas MPs | 0.9 | 1 = pure FLIP |

For a granular soil, a physically based choice of $\chi_{\text{off}}$ is the ratio of minimum to
current density: gasify when $e_p>e_{\max}$. This is the same idea as the existing
`$$MAXIMUM_POROSITY` input.

---

## 7. Implementation plan

### 7.1 New module and data

**Module `src/Gasification.FOR` (`ModGasification`):**

| Routine | Purpose |
| --- | --- |
| `InitialiseGas()` | set `V0Ref`, store `ESMstatev0Array`, set every MP to continuum |
| `UpdateGasState()` | once per step after the stress update: evaluate G1–G5 and R1–R2, update counters, switch states |
| `ReactivateParticle(p)` | §5.2 |
| `IsGasParticle(p)` | logical function used by the guards |
| `GasDiagnostics()` | number of gas MPs, switches per step, gas mass fraction |

**New particle fields** (`Particle.FOR`, with matching defaults and copies in `CopyParticle`):
`IsGas` (logical), `GasHoldCounter`, `RampCounter` (integer), `V0Ref` (real).
They are written to and read from the BRF restart (`WriteMPMData.FOR` and `ReadMPMData.FOR`), so the
state survives load steps and stages. `ESMstatev0Array` is persisted the same way.

### 7.2 Hooks in the existing code

| Step | File / location | Change |
| --- | --- | --- |
| Initialise | after MP initialisation (`MPMInit.FOR`) | `call InitialiseGas()` |
| State update | `MPMDYNConvPhase.FOR`, after `call MPMDYNGetSig()` (l.372) | `call UpdateGasState()` |
| Stress | `MPMDYNGetSig` (`MPMDYNStresses.FOR:274`, `:288`) | skip `StressSolid` for gas; apply $\lambda_r$ ramp |
| Stress mapping | `AssignStressesToParticles` (`MPMConvPhase.FOR:1025`) | skip gas MPs |
| Internal force | `MPMDYNBTSig` (`MPMDYNBTSig.FOR`) | skip gas MPs |
| Integration type | `CheckFillingOfElements` (`MPMConvPhase.FOR:952`) | gas → particle integration; fill counts continuum only |
| Damping | `MapDampingFromParticles` (`LagrangianPhase.FOR:2085`) | exclude gas MPs |
| Time step | `CalculateCriticalTimeStep` (`MPMDYNConvPhase.FOR:2882`) | skip wave speed for gas; add kinematic limit |
| Velocity | `UpdateParticleVelocityAndMapMomentum` (`MPMDYNConvPhase.FOR:946`) | optional $\alpha_g$ blending |
| Ghost | `BuildGhostFaces` | exclude elements holding only gas MPs |
| Output | `WriteVTKOutput.FOR` | MP scalar `gas_state` (0 continuum, 1 gas, 2 ramping) |

All guards test `CalParams%Gas%Apply` first, so results are bit-identical when the feature is off.

### 7.3 Inputs and problem type

**Kernel.** New derived type in `CalParams` (`ReadCalculationData.FOR`):

```fortran
type GasType
  logical :: Apply = .false.
  real(REAL_TYPE) :: StressTol = 0.01, PRef = 1.0
  real(REAL_TYPE) :: ChiOff = 0.60, ChiOn = 0.80, FillOff = 0.10, FillOn = 0.40
  integer(INTEGER_TYPE) :: NHold = 5, NRamp = 50
  real(REAL_TYPE) :: PICFactor = 0.9
end type
```

* Defaults are set in `InitialiseCalculationParameters` with `Apply = .false.`.
* Reader branch for `$$GASIFICATION` next to `$$APIC_FLIP_BLEND_FACTOR` (≈ l.1826), with `Assert`
  checks: $\varepsilon_\sigma>0$, $p_{\text{ref}}>0$, $0<\chi_{\text{off}}<\chi_{\text{on}}\le1$,
  $0\le\phi_{\text{off}}<\phi_{\text{on}}\le1$, $N_{\text{hold}},N_{\text{ramp}}\ge1$, $0\le\alpha_g\le1$.
* The kernel's CPS writer (`WriteCalculationParameters`, l.2887) re-emits the key for continuation files.
* Echo the active values to the log at stage start.

**CPS format:**

```
$$GASIFICATION
1  0.01  1.0  0.60  0.80  0.10  0.40  5  50  0.9
! flag eps_sigma p_ref[kPa] chi_off chi_on phi_off phi_on N_hold N_ramp alpha_g
```

**Per-material switch (GOM).** In `xml/Materials.xml`, add `allow gasification [No|Yes]` (default
Yes), hidden for the rigid body with the existing `hide_show_const_model` state proc. `createGOM.tcl`
writes `$$GASIFICATION_ALLOWED` after `$$MATERIAL_MODEL_SOLID`. `ReadMaterialData.FOR` reads it as an
optional key, and an absent key means allowed. `MatParams(I)%AllowGasification` feeds condition G5.

**Problem type** (`src/GiD_Problemtype/Anura3D_2025.gid`). The block sits in the same per-stage
container as ghost stabilisation, so it can be switched on and off per stage:

```
Calculation Data
└─ STABILISATION AND GASIFICATION
   └─ GASIFICATION   [do not apply gasification | apply gasification]
       ├─ gasification preset        [conservative | default | aggressive | user]
       ├─ stress tolerance [-]                 (user)
       ├─ reference pressure [kPa]
       ├─ density ratio to gasify [-]          (user)
       ├─ density ratio to reactivate [-]      (user)
       ├─ element fill to gasify [-]           (user)
       ├─ element fill to reactivate [-]       (user)
       ├─ hold steps [-]
       ├─ stiffness ramp steps [-]
       └─ gas FLIP factor [-]
```

`xml/CalculationData.xml`, inside the `STABILISATION` container defined in
`Ghost_Stabilisation_Anura3D.md` §6.3:

```xml
<value n="GASIFICATION" pn="GASIFICATION" icon="calculate.png"
       values="do not apply gasification,apply gasification"
       v="do not apply gasification" actualize_tree="1">
  <dependencies node="../value[@n='gas_preset']" att1="state" v1="normal" value="apply gasification"/>
  <dependencies node="../value[@n='gas_preset']" att1="state" v1="hidden" value="do not apply gasification"/>
  <!-- same pair for gas_p_ref, gas_hold_steps, gas_ramp_steps, gas_flip_factor -->
</value>
<value n="gas_preset" pn="gasification preset" values="conservative,default,aggressive,user" v="default"
       state="hidden" actualize_tree="1"
       help="conservative: gasify late, reactivate early; aggressive: gasify early; user: enter all values">
  <dependencies node="../value[@n='gas_stress_tol']" att1="state" v1="normal" value="user"/>
  <dependencies node="../value[@n='gas_stress_tol']" att1="state" v1="hidden" value="default"/>
  <!-- ... hidden for conservative/aggressive; same for the four threshold fields -->
</value>
<value n="gas_stress_tol"   pn="stress tolerance [-]"            v="0.01" state="hidden" help="|p|,q below tol*p_ref count as stress-free. 1e-3 to 0.1."/>
<value n="gas_p_ref"        pn="reference pressure [kPa]"        v="1.0"  state="hidden"/>
<value n="gas_chi_off"      pn="density ratio to gasify [-]"     v="0.60" state="hidden" help="rho/rho0 at or below this value gasifies. Must be below the reactivation ratio."/>
<value n="gas_chi_on"       pn="density ratio to reactivate [-]" v="0.80" state="hidden"/>
<value n="gas_fill_off"     pn="element fill to gasify [-]"      v="0.10" state="hidden"/>
<value n="gas_fill_on"      pn="element fill to reactivate [-]"  v="0.40" state="hidden"/>
<value n="gas_hold_steps"   pn="hold steps [-]"                  v="5"    state="hidden"/>
<value n="gas_ramp_steps"   pn="stiffness ramp steps [-]"        v="50"   state="hidden"/>
<value n="gas_flip_factor"  pn="gas FLIP factor [-]"             v="0.9"  state="hidden" help="1 = pure FLIP, lower values damp gas-particle noise."/>
```

`scripts/createCPS.tcl` (after the ghost block):

```tcl
# GASIFICATION
GiD_WriteCalculationFile puts {$$GASIFICATION}
set base {container[@n="Calculation_Data"]/container[@n="STABILISATION"]}
proc gasv {node base n} { return [$node selectNodes "string($base/value\[@n=\"$n\"\]/@v)"] }
set on     [gasv $stageNode $base GASIFICATION]
set preset [gasv $stageNode $base gas_preset]
switch $preset {
  conservative { set tol 0.005; set cOff 0.50; set cOn 0.85; set fOff 0.05; set fOn 0.50 }
  default      { set tol 0.01;  set cOff 0.60; set cOn 0.80; set fOff 0.10; set fOn 0.40 }
  aggressive   { set tol 0.05;  set cOff 0.75; set cOn 0.85; set fOff 0.20; set fOn 0.40 }
  user { set tol  [gasv $stageNode $base gas_stress_tol]
         set cOff [gasv $stageNode $base gas_chi_off];  set cOn [gasv $stageNode $base gas_chi_on]
         set fOff [gasv $stageNode $base gas_fill_off]; set fOn [gasv $stageNode $base gas_fill_on] }
}
set flag [expr {$on eq "apply gasification" ? 1 : 0}]
GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s %s %s %s %s" $flag $tol \
    [gasv $stageNode $base gas_p_ref] $cOff $cOn $fOff $fOn \
    [gasv $stageNode $base gas_hold_steps] [gasv $stageNode $base gas_ramp_steps] \
    [gasv $stageNode $base gas_flip_factor]]
```

The preset numbers above are starting proposals, to be fixed by the tuning runs in §8.

---

## 8. Diagnostics and verification

**Per-step log line:** number of gas MPs, gasified and reactivated counts this step, gas mass
fraction, number of elements holding only gas MPs, the controlling time step and which limit set it
(wave speed or kinematic), and the total-momentum change.

**Benchmarks:**

1. **Feature off.** Results must be bit-identical to the current build.
2. **Ballistic block.** A small stress-free block launched upward should follow
   $\mathbf x(t)=\mathbf x_0+\mathbf v_0t+\tfrac12\mathbf g t^2$ exactly once gasified, with momentum
   change equal to $m\mathbf g\,\Delta t$ per step.
3. **Drop and redeposit.** A loose column falls onto a rigid base. Check that MPs reactivate, the
   deposit reaches a stress state that grows with depth after the ramp, and no element inverts.
4. **Granular column collapse.** Compare runout, number of active elements and run time with the
   feature on and off.
5. **`Free_fall_wall_BH` cone model.** Fewer active elements, no "outside the mesh" errors, a time
   step no longer controlled by ejecta elements, and a cone reaction force within a few percent of the
   reference where the reference is stable.
6. **Threshold sweep** over $\chi_{\text{off}}$, $\phi_{\text{off}}$ and $\varepsilon_\sigma$, which
   also confirms there is no chattering, i.e. few MPs switching state more than once.

---

## 9. As implemented

| Item | Implementation |
| --- | --- |
| Module | `src/Gasification.FOR` (`ModGasification`). State is kept in module arrays indexed by MP (`GasState`, `GasHold`, `GasRamp`, `GasV0`), not in `Particle`. |
| Input | `$$GASIFICATION` in the CPS: `flag tol pRef chiOff chiOn fillOff fillOn nHold nRamp gasFLIP`. Off by default. It is re-written to the continuation CPS files. |
| Element fill | $\phi_K=\sum_{p\in K}V_{p,0}/\lvert K\rvert$, i.e. **material content** from the reference volumes, gas and continuum MPs together. The current volume cannot be used: on the cone model, flying MPs swell 6–20×, so a lone MP showed a fill of 0.5–1.4 instead of about 1/12. |
| Continuum → gas | Element not blocked (no rigid-body or other-entity MPs, no contact interface nodes) and, for $N_{\text{hold}}$ steps, any of the following. **(a)** Stress-free (G1, G2), loosened or sparse (G3), and $\dot\varepsilon_v\ge0$. **(b)** Isolated: $\phi_K\le\phi_{\text{off}}$ in the element and in all of its edge neighbours. **(c)** Loose: $\chi_p\le\chi_{\text{off}}$ **and** $\phi_K\le\phi_{\text{off}}$. Branches (b) and (c) ignore the MP stress, because a lone MP in a sparse element is traction free and the stress it carries is spurious; on the cone model such MPs carried 5–25 kPa. |
| Gas → continuum | R1 or R2 for $N_{\text{hold}}$ steps. Re-initialisation: stress set to 0, $V_{p,0}\leftarrow V_p$, ESM state variables restored from a copy taken at the first call, then a strain-increment ramp over $N_{\text{ramp}}$ steps (applied only in MP-integrated elements). The porosity is left to `DynUpdateParticlePorosity`. |
| Gas MP | `MPMDYNGetSig` sets $\boldsymbol\sigma=0$ and skips `StressSolid`. Elements holding gas MPs are forced to MP integration after `CheckFillingOfElements`. Gas MPs are excluded from local damping, from the wave-speed time step and from the ghost-element fill and density. They use their own FLIP fraction in `UpdateParticleVelocityAndMapMomentum`. |
| Time step | Kinematic bound $\Delta t\le h_K/\lvert\mathbf v_p\rvert$ for gas MPs (times the Courant number). |
| Restart | BRF block `$$GASSTATE$$` (4 reals per MP), written only when gasification is on and read by `ReadGasStateFromFile` at the start of the next stage. |
| Output | VTK MP scalar `gas_state` (0 continuum, 1 gas, 2 ramping). Log line once per load step: `Gas: points, gasified, reactivated, ramping, candidates isolated/loose/stress-free/all, mass fraction, gas-only elements, kinematic dt`. |
| GUI | `GASIFICATION` toggle, preset and fields in the `STABILISATION` container (`CalculationData.xml`), written by `createCPS.tcl`. Older models write it as off. |
| Not done yet | Per-material switch (§7.2, `$$GASIFICATION_ALLOWED`). Rigid-body MPs are always excluded. 3D TET4 is untested. Benchmarks 2–4 and 6 of §8 have not been run. |

**Cone model (`Free_fall_wall_BH`, 25 load steps of 0.005 s, ghost 0.2/0.2/0.6, default gas preset).**

| Load step | Ghost only: max soil speed | Ghost + gas: max soil speed | Gas MPs | Cone tip $y$ (ref / gas) |
| --- | --- | --- | --- | --- |
| 10 | 18.4 m/s | 18.2 m/s | 74 | −0.1159 / −0.1157 |
| 13 | 17.5 m/s | 17.2 m/s | 89 | −0.1241 / −0.1237 |
| 14 | **262 m/s** | 16.6 m/s | 99 | −0.1269 / −0.1265 |
| 15 | **8378 m/s**, then "Maximum number of recursions reached" | 16.2 m/s | 120 | – / −0.1293 |
| 25 | – | 17.6 m/s | 294 | – / −0.1372 |

Gas MPs reach at most 1–4 m/s. The gas mass fraction is below 1 %, and some MPs are reactivated in most load steps. The time step stays at about 1e-4 s. With gasification on but no MP gasified, the log is identical to the ghost-only run.

