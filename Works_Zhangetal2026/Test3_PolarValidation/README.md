# Polar-axis analytical validation

Four sources (normal, small, 10-m thin and full-longitude zonal band) are evaluated at 2,001 exact-North-Pole observations. Eight selected quantities are compared with an analytical reference.

`PlotTest3_AllFigures.m` can plot archived `results/polar_validation.mat` from the data attachment without recalculating the analytical solution.

For a fresh run, build CMake target `polar_cpp_validation` with research ON, create `results/`, and run it with this directory as working directory. Configure `TLI_POLAR_REFERENCE_DIR` to the external `Tesseroid_AnalyticalPolar_vpa.m` directory, then run `Test3_PolarValidation.m` in MATLAB with Symbolic Math Toolbox. That analytical-reference function was not in the supplied source directory and is not bundled.

The exact-pole numerical branch is distinct from ordinary off-axis line quadrature. See [method](../../docs/METHOD.md) and [reproduction](../../docs/REPRODUCIBILITY.md).
