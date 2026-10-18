# Laplace–Poisson validation

The data attachment supplies the three `TFG_<height>.ForPar` inputs and `ResultsPan_<height>.dat` outputs at +1, -31 and -15 km, each with 160,801 observations. All use legacy scaled units. The driver `Cpp/Test4_Work1_GravityEstimate.cpp` contains fixed filenames: select the intended input and output before each new run. Build it with root CMake target `trace_validation` (research ON).

`AnalyzeLaplaceError.m` and the `Figure5/` scripts analyze archived trace residuals; the unified `Analysis/RunValidationAnalysis.m` computes the manuscript statistics. Optional `CompareTracePan.m` additionally requires external historical comparison outputs configured by `TLI_TRACE_REFERENCE_DIR`.

Trace residuals test differential consistency, not independent accuracy of all tensor components. Near material boundaries and for strict interior high derivatives, independent componentwise evidence remains a separate requirement. See [reproduction](../../docs/REPRODUCIBILITY.md).
