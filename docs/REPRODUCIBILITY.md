# Reproducing experiments

## Distribution layers

The source ZIP is suitable for committing to GitHub. It contains production code, scripts, small CSV summaries, release checks and documentation. The separately supplied data ZIP is a release attachment: extract it into the repository root. Its paths begin with `Works_Zhangetal2026/` or `Tests/`, and `docs/DATA_MANIFEST.csv` lists its original file hashes and sizes.

Neither ZIP contains the manuscript drafts, PDF/reference copies, generated bitmap figures, precompiled executables/MEX, editor backups, global source-model `.ForPar` files, or external third-party assets. Original-file inventory and distribution decisions are recorded in [SOURCE_INVENTORY.csv](SOURCE_INVENTORY.csv).

Source is MIT licensed. Archived numerical result arrays in the data attachment are supplied as accompanying research outputs; externally sourced ETOPO/VIM model content and map assets are excluded. A downstream data repository should record a separate data license chosen by the authors; no independent data-license grant is invented here.

## Fast release checks

Run the CMake/CTest commands in [BUILD.md](BUILD.md). These tests use no external datasets. They compare all 20 ordinary-TLI components with independently integrated Cartesian volume kernels at two separated nonsymmetric observations, check the reference's 20-to-32-node convergence, check the interior Poisson trace and rejected geometry, and compare a constant global TLI-MMSI shell with closed-form results in its three radial regimes. The original regional-grid finite-value smoke test is also retained.

The direct volume reference is intentionally used only at separated exterior points. It would miss distributional contact terms if applied naively inside the source.

## Manuscript experiments

Run scripts from their own experiment directories unless they determine their location from `mfilename`. Preserve archive data before recomputing: the historical drivers use fixed output filenames. Full numerical runs and high-precision quadrature can take hours or require a cluster.

| Folder | Plot / analyze archived data | Recompute |
| --- | --- | --- |
| `Test1_GaugeDimensionalReduction` | `figures/plot_Figure2_DimensionalReduction.m` reads the two `summary_current.csv` files | Build MEX; run `run_plus100km`, `run_plus400km`; requires Symbolic Math Toolbox; checkpoints resume existing runs |
| `Test2_SphericalShellValidation` | `Figures/plot_plus1km`, `plot_minus31km`, `plot_minus15km` read the three archived `.dat` files | Build CMake `shell_validation` with MPI/research ON; run from this folder, e.g. `mpirun -np 4 /path/to/shell_validation 60` for sampling; stride 1 is full grid |
| `Test3_PolarValidation` | `PlotTest3_AllFigures.m` reads `results/polar_validation.mat` | Build/run `polar_cpp_validation` from this folder with `results/` present, then run `Test3_PolarValidation.m` with the external analytical function configured |
| `Test4_LaplacePoissonValidation` | `Figure5/Test4_Plot*.m`, `AnalyzeLaplaceError.m`; `.dat` files in data attachment | Build `trace_validation`; inspect/edit its fixed parameter/result filenames per height; three `.ForPar` inputs are included in the data attachment |
| `Analysis` | `RunValidationAnalysis.m` reads Test1 summaries, polar MAT and Test2/Test4 results | Regenerates CSV statistics and TeX tables from the supplied archived arrays; does not rerun forward modelling |
| `RunTimeTest1` | `Analyze_RuntimeTest1.m`, `Figures/Figure_RuntimeTest_1km.m` read the full timing CSV | Build `Test_RuntimeTest` with MPI/research ON; full run has 160,801 observations and repetitions |
| `RunTimeTest2` | Plotting requires archived result files and external map resources | Generate `MSIG.ForPar` from external resampled ETOPO arrays, then run appropriate mass-shell MPI driver |
| `RunTimeTest3` | Plotting requires `ResultM.dat` and external map resources | Generate `MSIM.ForPar` from external VIM data, then run appropriate magnetization-shell MPI driver |

MATLAB R2025b was used in the source project. Symbolic Math Toolbox is needed for VPA references; Mapping Toolbox and external plotting helpers are needed by the global map scripts. MATLAB's own OpenMP MEX wrappers do not need Parallel Computing Toolbox.

## External inputs

The following dependencies were absent from the requested source directory. Their paths are now environment variables, with explicit missing-input messages:

| Variable | Expected contents | Needed by |
| --- | --- | --- |
| `TLI_ETOPO_DIR` | `ETOPO_<spacing>Degree.dat` tables in the exact resampled layout expected by the scripts; height column is converted from km to m | RunTimeTest2 input generator; coordinate generator writes here |
| `TLI_VIM_FILE` | Path to `7-VIMModel.dat`, in the source project's VIM column convention | RunTimeTest3 input generator |
| `TLI_POLAR_REFERENCE_DIR` | Directory containing `Tesseroid_AnalyticalPolar_vpa.m` and its dependencies | Recomputing polar analytical reference |
| `TLI_TRACE_REFERENCE_DIR` | External historical comparison `ResultsPan_<height>.dat` files | Optional `CompareTracePan.m` only |
| `TLI_GEOCONTOUR_DIR` | `SAGP_Contour` helpers, shapefiles and MAT resources named in the map scripts | Recreating the original global map style |

Example configuration in MATLAB:

```matlab
setenv('TLI_VIM_FILE','/path/to/7-VIMModel.dat')
setenv('TLI_POLAR_REFERENCE_DIR','/path/to/analytical-reference-functions')
```

The release does not invent download URLs, provenance, original dataset version numbers or redistribution permissions for these local inputs. To make global benchmark input generation fully self-contained, the authors should archive the resampling scripts and precise external model references. The provided raw benchmark outputs still allow result analysis, with the limitations stated below.

## Version correspondence

Several old README files described an earlier far-field GL6 line path, whereas the current production `point_estimate` uses 2D far-field integration. Older runtime notes also described a different source discretization. The manuscript itself records this provenance gap. Release documentation now describes the current source, while copied numerical tables and archived binary arrays remain historical evidence. Their original SHA-256 values are retained. They have **not** been relabelled as freshly recomputed results of this publication candidate.

Before attaching scientific accuracy claims to a final tag, recompute the relevant full experiments at that tag and record compiler, flags, hardware, effective tolerances, input hashes, convergence/budget diagnostics, output hashes and timings. In particular, constant-shell symmetry and small trace residuals do not replace componentwise checks for nonuniform shells or near-boundary single-body fields.
