# Release validation

Executed locally on **2026-10-18**, Windows x64. These results apply to the publication candidate, not to every archived experiment.

## Environment

- MSVC 19.51.36260.0, Visual Studio 2026; C++17, Release configuration.
- CMake 4.3.1-msvc1, CTest; OpenMP 2.0 as detected by CMake.
- Intel MPI 2021.17 (MPI interface version 4.1), two local ranks for driver checks.
- MATLAB R2025b; MEX rebuilt with the supplied Intel oneAPI build script and OpenMP.
- The installed MATLAB did not detect a supported compiler through standard `mex`; the portable MATLAB-selected-compiler builder therefore remains untested locally. The Intel fallback was compiled and executed successfully.

## Newly executed numerical checks

| Check | Result |
| --- | --- |
| All 20 ordinary-TLI components vs independent Cartesian tensor-product volume integrals at two nonsymmetric exterior observations | Passed; maximum componentwise relative discrepancies `5.8674e-13` and `3.80255e-13` |
| Independent reference refinement from 20 to 32 Gauss nodes per dimension | Passed the release reference-stability gate |
| Strict-interior Poisson trace, zero-density model, invalid radii and undefined lower-face rejection | Passed |
| Original finite-value TLI-MMSI smoke test | Passed |
| Constant global shell vs closed forms outside, in cavity and in material | Passed for all 20 components at three sampled radii |
| CTest suite | 3/3 passed; observed total 22.03 s on this machine |
| Ordinary gravity and magnetic serial vs OpenMP vs two-rank MPI batch outputs | Identical in the exercised cases |
| Shell gravity all-components vs Vz-only, two ranks | Identical for selected output |
| Shell magnetic all-components vs vector-only, two ranks | Maximum relative difference `5.65e-15` for the selected vector |
| MATLAB serial/OpenMP gravity and magnetic, superposition, 20-output canonical/public mapping and km-to-SI conversion vs separately built C++ SI example | Passed |
| MATLAB negative-tolerance validation error identifier | Passed after the gateway diagnostic fix |
| MATLAB small demonstration | Passed, finite gravity/magnetic outputs |

The independent-volume acceptance gate is a relative `2e-7`; observed discrepancies above are much smaller. The constant-shell test uses a derivative-order-scaled absolute gate of `2e-7` and a coarse 10-degree source grid. These finite sets are regression evidence, not an exhaustive proof near every boundary or for nonuniform shells.

All six production MPI targets plus `shell_validation` and `Test_RuntimeTest` compiled. All default/ordinary OpenMP targets and the optional polar, trace and kernel-shape research targets compiled. Full research executables were not run as part of those build checks.

## Archived-data analysis

Using isolated working copies of the supplied arrays, MATLAB successfully reran:

- `Analysis/RunValidationAnalysis.m` for Tests No. 1–11;
- `RunTimeTest1/Analyze_RuntimeTest1.m` for the 160,801-point timing table.

This verifies that the retained arrays and audit scripts can be read and analyzed together. It does not recompute those forward calculations or establish their exact historical source version. The source archive retains the original compact summary tables rather than overwriting them with the audit copy.

Five regenerated validation/runtime CSV tables matched the original tables exactly. Both ZIP archives passed CRC checks; extraction followed by the shipped verifier checked 191 source/data files against SHA-256 manifests. Local Markdown link targets were checked. All production header hashes match the supplied originals. GitHub workflow and citation metadata were inspected, but a CFF schema validator was not installed.

## Limits of this validation

Linux/macOS and GitHub-hosted CI have not executed in this local session. The workflow is provided for that follow-up. Full VPA, 64,800-point shell, complete exact-pole profiles, 160,801-point forward trace experiments and global ETOPO/VIM forward benchmarks were not rerun. No new universal speedup or machine-precision claim is made for the hybrid algorithm. External data dependencies and the historical source-version gap are documented in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
