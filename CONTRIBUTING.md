# Contributing

Please include geometry, units, component names, tolerances, compiler/build flags, convergence diagnostics and a minimal reproducer when reporting numerical issues. State whether the issue concerns ordinary TLI or the hybrid shell model.

Run the CMake/CTest checks before submitting changes. Numerical-kernel changes need an independent reference with a clear error scale, not only trace/symmetry checks. Preserve the distinction between the ordinary and shell kernel variants. Do not commit local executables, MEX binaries, generated maps or large result arrays.
