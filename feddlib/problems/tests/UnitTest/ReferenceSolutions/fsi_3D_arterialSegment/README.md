# Frozen arterial outlet restart references

These bundles were generated from uninterrupted four-rank simulations of the
supplied 0.5 mm artery, using the settings in
`RestartTests/fsi_3D_arterialSegment`. They are numerical regression baselines;
normal tests only read them and never regenerate them.

Each model directory contains complete restart state, not just final solution
vectors: coupled/component manifests, primary fields, ALE geometry, fluid
mass/solution history, Newmark displacement/velocity/acceleration, multiplier
fields and `FSIOutletState_<time>.xml` snapshots. Metadata uses checkpoint format
version 2, including actual clock counters, dt/dt_prev and history timestamps.

| Model | Stored snapshots | Normal restart comparisons |
|---|---|---|
| Resistance | Every step from 0 to 0.015 | 0.005 -> 0.010 and 0.010 -> 0.015; save-all chain 0.005 -> 0.010, then 0.0075 -> 0.0125 |
| Absorbing | 0.005, 0.010, 0.015 with prior history | 0.005 -> 0.010 and 0.010 -> 0.015 |
| Absorbing Paper | 0.005, 0.010, 0.015 with prior history | 0.005 -> 0.010 and 0.010 -> 0.015 |

The reference producer stops at 0.0175 so the existing FSI start-of-step export
writes the 0.015 snapshot. Normal tests advance exactly two timesteps per restart
with dt=0.0025. The transition time 0.00625 lies within the first two-step case;
the later case checks restoration of an already captured transition area.

`Generation.txt` and `GenerationInputs` record the producer settings and
command-line overrides. The nonlinear and linear tolerances are both 1e-10;
the field comparison retains `error_l2 <= 1e-15 + 1e-12 * norm(reference)` and
the scalar history bound remains `1e-14 + 1e-12 * abs(reference)`.

Regenerate deliberately with the script documented in
[`fsi_3D_arterialSegment/README.md`](../../../RestartTests/fsi_3D_arterialSegment/README.md).
Review numerical changes before replacing these source bundles, then rebuild.
CTest mutates only isolated copies for missing-state/version compatibility tests.
