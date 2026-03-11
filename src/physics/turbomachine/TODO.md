# Turbomachine TODO

## Near-term modeling work

- Update the blade aero model to handle blade stall causally, not just diagnostically.
  - Current state: `stall_margin` and `stall` are reported, but a stalled row still marches with the same basic turning/loss closure.
  - Target: make stall affect turning, loss, or row validity so off-design/stalled points stop looking artificially well-behaved.

- Add explicit post-stall policy to the axial streamtube model.
  - Options to evaluate:
    - hard invalidation once stall limit is exceeded
    - soft degradation of turning and increased loss beyond stall onset
  - The first pass can be hard invalidation if we want stall to define the feasible operating envelope.

- Tighten the physical reconstruction used by axial operating-point diagnostics.
  - Current station `mdot_error` shows that EOS-based reconstruction is not perfectly self-consistent with the streamtube model’s internal nondimensional closure.
  - Add relative mass-flow error reporting and decide whether reconstruction or solver closure should be adjusted.

- Improve the demo turbine design using the new diagnostics.
  - Tune metal angles and deviation parameters so the turbine demo has:
    - expansion over a useful operating region
    - non-stalled nominal operating point
    - sane efficiency and loss distribution row-by-row

- Evaluate dimensionalizing the axial streamtube model.
  - Motivation:
    - reduce dependence on ideal-gas-style reference scaling
    - improve EOS consistency between the solver and the physical reconstruction
    - make diagnostics and continuity interpretation more direct
  - Important constraint:
    - dimensionalization will not, by itself, remove the implicit outlet-velocity solve or multi-branch behavior
    - it should be treated as a later architectural cleanup, not as the first fix for current turbine-map artifacts

## Diagnostics and analysis

- Add a diagnostic summary of first-failing row and failure mode.
  - Examples:
    - first stalled row
    - first choked row
    - first invalid row

- Add a direct comparison utility between map-based operating-point results and replayed streamtube results over a sweep.
  - Goal: quantify where the tabulated map departs from the source axial model.

- Consider CSV export for `diagnose_axial_operating_point.jl`.
  - Separate station and row tables would make it easier to inspect and post-process.

## Meanline / map generation

- Decide whether feasible-map boundaries should be limited by stall once causal stall handling is implemented.

- Revisit map tabulation defaults for turbine demos after the turbine axial model is retuned.
  - Current plots are useful, but the demo turbine is still not a strong reference case.

## Structural cleanup

- Remove remaining compatibility aliases that still use `meanline` naming now that the common API is axial-model based.

- Decide whether compressor and turbine residual equations can be merged now that the canonical map is `TabulatedPerformanceMap`.
