# Manuscript experiments

These folders retain experiment sources and small summaries for the Zhang–Du tesseroid line-integration manuscript. Large numerical arrays are distributed in a separate data ZIP.

See [the reproduction guide](../docs/REPRODUCIBILITY.md) for dependencies, working directories, rebuild commands and the distinction between archived results and newly executed checks. See [the method guide](../docs/METHOD.md) for the current hybrid algorithm.

| Folder | Experiment |
| --- | --- |
| Test1_GaugeDimensionalReduction | Independent high-precision 3D / 2D / 1D comparisons for 20 components at two exterior points |
| Test2_SphericalShellValidation | Complete homogeneous shell, three radial regimes |
| Test3_PolarValidation | Four source geometries on the exact polar axis |
| Test4_LaplacePoissonValidation | Trace and differentiated trace residuals |
| RunTimeTest1 | Single-tesseroid MPI timing |
| RunTimeTest2 | Global topographic-shell timing |
| RunTimeTest3 | Global magnetization-shell timing |
| Analysis | Statistics recomputed from archived validation arrays |

The snapshot lacks an exact source-version link for some archived runs. Those arrays remain historical research data. Global source models, resampling provenance and several external plotting/reference dependencies are not contained in the original requested directory. The package does not claim all experiments are self-contained.
