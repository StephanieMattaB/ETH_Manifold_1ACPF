# Legacy / superseded scripts

These scripts were an earlier, exploratory pass at the DSO pipeline. They are
kept for reference and git history but are **superseded** by the consolidated
pipeline at the repository root (`run_ieee13_pipeline.m` and its helpers:
`bd_linear_solve.m`, `bd_tangent_general.m`, `ieee13_der.m`,
`ieee13_congest.m`, `ieee13_capacitors.m`, `voltage_report.m`,
`compare_voltages.m`, `derive_shared_constraint.m`, `run_stage.m`).

Known issues in this legacy code, fixed in the current pipeline:
- `run_ieee13_nonlinear_ac.m` hardcoded `stress_factor = 1.3`, silently
  inconsistent with `example_ieee13_congested_variant.m`'s `load_scale = 1.8`
  (both are unified on `lambda = 1.8` in the current pipeline).
- There was no separate active-system stage `S'` (inverter-based DER at
  645/611/652/671): congestion was applied directly to the original
  benchmark `S`, skipping the `S -> S' -> S''` structure.
- `build_shared_constraint` was duplicated (once embedded in
  `example_ieee13_congested_variant.m`, once in `relinearize_ieee13_congested.m`,
  plus a third generic/untied copy under `../constraint/`), and was invoked at
  more than one stage rather than only at the final `S''`.

`example_ieee13_congested_variant.m` is the modified version of the paper's
original `example_ieee13.m` example that this branch had been iterating on;
the root `example_ieee13.m` has been restored to the original, unmodified
paper example.

`../data/network_to_dynect.mat` (moved here) and
`nonlinear_ac_validation_results.mat` are stale outputs of this legacy code.
