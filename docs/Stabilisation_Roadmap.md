# Roadmap: ghost stabilisation and gasification in Anura3D

**Status:** plan for review. No code has been written yet.
**Related documents:** `Ghost_Stabilisation_Anura3D.md`, `Gasification_Anura3D.md`.
**Scope:** explicit dynamic MPM-MIXED, dry single-phase, TRI3 (2D plane strain / axisymmetric) first,
TET4 after.

---

## 1. What each feature fixes

| Problem seen in the impact runs | Cause | Feature |
| --- | --- | --- |
| Velocity spikes and stress oscillations at free surfaces; MPs thrown "outside the mesh" | nodes with tiny lumped mass in cut (partially filled) elements | **Ghost stabilisation** |
| Many nearly empty elements from ejecta; time step and run time driven by them | sparse, stress-free MPs treated as a stiff continuum | **Gasification** |

The two features are complementary. Gasification removes stiffness from MPs that no longer behave as
a continuum. Ghost stabilisation makes the edge of the remaining continuum well behaved.

---

## 2. Recommended order

### Phase 0: diagnostics (no physics change)

Add a per-step log line and VTK fields for:
* the number of active, particle-integrated and nearly empty elements ($\phi_K<0.1$);
* the minimum nodal mass ratio $\min_a m_a/\bar m$;
* the controlling time step and which element or node sets it;
* total linear momentum and its change per step.

These numbers give the baseline that both features are judged against. Run them on the
`Free_fall_wall_BH` cone model and keep the output as the reference.

### Phase 1: ghost stabilisation

| Step | Content | Exit test |
| --- | --- | --- |
| 1a | GUI toggle and presets, CPS key `$$GHOST_STABILISATION`, kernel reader (default off), face detection and `ghost_element` output only | feature off is bit-identical; ghost faces look right in ParaView |
| 1b | ghost force (stiffness / viscosity term) | momentum change exactly zero; elastic-column test improves |
| 1c | ghost mass: diagonal option first, then Jacobi | minimum nodal mass ratio raised; momentum error below $10^{-10}$ with Jacobi |
| 1d | ghost time-step bounds and final presets from a parameter sweep | stable time step independent of fill; presets fixed |

### Phase 2: gasification

| Step | Content | Exit test |
| --- | --- | --- |
| 2a | GUI toggle, per-material switch, CPS key `$$GASIFICATION`, reader (default off); **classification and `gas_state` output only**, with no change in behaviour | thresholds tuned on the cone model: the MPs flagged as gas are the ejecta |
| 2b | gas behaviour: zero stress, no internal force, damping exclusion, integration type, kinematic time step | ballistic-block test exact; momentum conserved |
| 2c | reactivation: density reset, state-variable reset, stiffness ramp | drop-and-redeposit test stable; no chattering |
| 2d | BRF restart persistence, continuation CPS re-emission, final presets | multi-stage run keeps the gas state across stages |

### Phase 3: coupling and tuning

* Ghost faces exclude elements that hold only gas MPs. Mixed elements use their continuum fill.
* Run a joint sweep on the cone model (presets × presets), then fix the GUI defaults.
* Write a short user note for the problem-type help.

In each phase the GUI and CPS plumbing comes first, so every later step can be tuned from the
problem type without recompiling.

---

## 3. Why ghost stabilisation comes first

1. **Reactivation needs it.** When gas MPs land and switch back, they create freshly filled, partially
   filled elements at the deposit surface. Without ghost stabilisation, these are exactly the unstable
   cut elements, and gasification would trade one instability for another.
2. **It can be verified on its own.** It is deterministic, has published benchmarks (Coombs 2023) and
   exact conservation properties (zero row sums), so errors are easy to find.
3. **It builds shared infrastructure.** Fill classification, face lists, diagnostics and the
   `STABILISATION` GUI container are reused by gasification.
4. **Gasification thresholds only make sense on a stable continuum.** If the continuum edge is still
   noisy, the stress-free and density tests in G1–G3 trigger on noise, not on real ejecta.

Gasification step 2a (classification only) can be developed **in parallel** with phase 1, because it
changes no behaviour. Its output helps choose the thresholds early.

---

## 4. Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Ghost penalty too strong: over-stiff edges, small time step | viscous bound in the time step; presets from a sweep; log the controlling limit |
| Diagonal ghost mass breaks momentum conservation | use the Jacobi option as default; log the momentum error |
| Gas MPs chatter between states | hysteresis ($\chi_{\text{off}}<\chi_{\text{on}}$, $\phi_{\text{off}}<\phi_{\text{on}}$) plus the hold counter |
| Shock on reactivation | stress-free reset and stiffness ramp |
| Gas state lost between load steps | BRF persistence and CPS re-emission (step 2d) |
| Interaction with the Bardenhagen contact | ghost terms never touch the rigid entity; gas MPs stay soil MPs, so contact still applies |
| Old models change behaviour | both features default off; bit-identical test in every phase |

---

## 5. Benchmarks (shared)

1. Elastic column under self-weight, offset on the grid (ghost).
2. Vibrating elastic bar with a partially filled end element (ghost: energy, momentum, time step).
3. Ballistic block (gasification: exact trajectory).
4. Drop and redeposit (gasification: reactivation).
5. Granular column collapse (both: runout, active elements, run time).
6. `Free_fall_wall_BH` cone impact (both: stability, active elements, cone reaction force).
