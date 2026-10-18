# Dimensional-reduction validation

At 100 km and 400 km above the upper source face, evaluate all 20 components with separately evaluated VPA volume, radially reduced surface and gauge-line calculations, then compare with the production C++ MEX implementation. Geometry is nonsymmetric.

Build the MEX and install Symbolic Math Toolbox. Run `run_plus100km.m` and `run_plus400km.m`. Existing checkpoints resume; use a separate checkout or preserve the archived results for a fresh run. Working precisions, convergence gates and integration requests are defined in `test1_settings.m` and `test1_run_case.m`.

`figures/plot_Figure2_DimensionalReduction.m` reads the compact `results/plus100km/summary_current.csv` and `results/plus400km/summary_current.csv` files included in the source archive. Checkpoints and detailed history are in the data attachment. A smaller cross-path discrepancy than the quadrature request is an observed comparison, not a certified bound.

See [reproduction details](../../docs/REPRODUCIBILITY.md).
