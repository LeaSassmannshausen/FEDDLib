# Frozen arterial outlet restart references

These bundles were generated from uninterrupted four-rank simulations of the
supplied 0.5 mm artery, using the settings in
`restartTests/fsi_3D_arterialSegment`. They are numerical regression baselines;
normal tests only read them and never regenerate them.

Each model directory contains complete restart state, not just final solution
vectors: coupled/component manifests, primary fields, ALE geometry, fluid
mass/solution history, Newmark displacement/velocity/acceleration, multiplier
fields and `FSIOutletState_<time>.xml` snapshots. Metadata uses checkpoint format
version 2, including actual clock counters, dt/dt_prev and history timestamps.

| Model | Stored snapshots | Normal restart comparisons |
|---|---|---|
| Resistance | Every step from 0 to 0.015 | One restart, 0.005 -> 0.010 |
| Absorbing | 0.005, 0.010, 0.015 with prior history | One restart, 0.010 -> 0.015 |
| Absorbing Paper | 0.005, 0.010, 0.015 with prior history | One unsteady-phase restart, 0.010 -> 0.015 |

The existing Resistance bundle was generated with save-all output and contains
additional intermediate snapshots. The current restart test reads its frozen
checkpoint at 0.005; future reference generation uses selected checkpoints for
all three models.

The reference producer stops at 0.0175 so the existing FSI start-of-step export
writes the 0.015 snapshot. Normal tests advance exactly two timesteps per restart
with dt=0.0025. Only Absorbing Paper has an outlet transition, at 0.00625;
its restart checks restoration of an already captured transition area.

`Generation.txt` and `GenerationInputs` record the producer settings and
command-line overrides. The bundles were regenerated on 2026-10-21 using the
current inputs: nonlinear tolerance 1e-6, coupled linear tolerance 1e-6,
and geometry linear tolerance 1e-10;
the field comparison retains `error_l2 <= 1e-15 + 1e-12 * norm(reference)` and
the scalar history bound remains `1e-14 + 1e-12 * abs(reference)`.

Regenerate deliberately with the script documented in
[`fsi_3D_arterialSegment/README.md`](../../../restartTests/fsi_3D_arterialSegment/README.md).
Review numerical changes before replacing these source bundles, then rebuild.
CTest mutates only isolated copies for missing-state/version compatibility tests.
