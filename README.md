# Tesseroid Line Integration

C++17 and MATLAB code for gravity and magnetic forward modelling in spherical geometry. Ordinary **TLI** evaluates the potential and Cartesian derivatives through third order for a homogeneous tesseroid. **TLI-MMSI** combines near-region line integration with far-region surface integration for gridded mass and magnetization shells.

[中文说明](README.zh-CN.md) · [Build](docs/BUILD.md) · [API and units](docs/API.md) · [Method](docs/METHOD.md) · [Reproduction](docs/REPRODUCIBILITY.md) · [Validation](docs/VALIDATION.md)

This is the `0.1.0-rc.1` publication candidate, prepared from the authors' local research snapshot on 2026-10-18. The accompanying manuscript is **“An exact line-integral formulation for evaluating tesseroid potentials and their derivatives through third order”**, Pan Zhang and Jinsong Du. The public repository is https://github.com/ZhangPan-666/Tesseroid-LineIntegration. A publication DOI has not yet been assigned.

## Quick start

Requires a C++17 compiler and CMake 3.18 or newer. The default build has no MATLAB, MPI, OpenMP, Eigen, or external-data requirement.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
```

Run `build/tli_example` on a single-configuration generator, or `build/Release/tli_example.exe` with Visual Studio. It prints all 20 components, estimated errors, and convergence flags for a small SI-unit gravity example. See [examples/gravity.cpp](examples/gravity.cpp).

For direct inclusion, add `C++Codes/` to your include path and include `_TLI_estimate.h`. With CMake, use `add_subdirectory(...)` and link `Tesseroid::tli`. The separate shell target is `Tesseroid::mmsi`.

## MATLAB

```matlab
addpath('Functions')
build_Tesseroid_LineIntegration_mex
model = [34.8 35.2 119.7 120.3 6341.2 6371.2 2670];
points = [35.11 120.08 6771.2];
[V,Vx,Vy,Vz] = Tesseroid_GravityEstimate(model,points,1e-12,1e-9);
```

The MATLAB/Work compatibility interfaces use **degrees and G = 1**. The example above uses kilometres and returns scaled quantities. Multiply by the appropriate physical constant and convert length units as described in [API.md](docs/API.md). The direct C++ API uses **radians** and defaults to SI `G = 6.67430e-11`. These conventions must be kept distinct.

The `P` wrappers request OpenMP inside the MEX file. They do not use MATLAB `parfor` or require Parallel Computing Toolbox. Serial builds can be requested with `build_Tesseroid_LineIntegration_mex(false)`.

## Contents

| Directory | Purpose |
| --- | --- |
| `C++Codes/` | Ordinary TLI headers, serial/OpenMP/MPI drivers and binary-I/O MATLAB scripts |
| `TLI_MMSI/` | Hybrid gridded-shell headers and four pure-MPI drivers |
| `Functions/` | MEX gateway, serial and OpenMP MATLAB wrappers, build scripts |
| `examples/` | Small gravity and MATLAB examples requiring no external data |
| `Tests/` | Fast release checks and optional kernel-shape research programs |
| `Works_Zhangetal2026/` | Manuscript experiment sources, plotting scripts and compact summaries |
| `C++Codes_Cluster/` | Retained standalone ordinary-TLI MPI layout; CMake uses `C++Codes/` |
| `docs/` | Method, conventions, reproduction, provenance and release notes |

Large numerical results are supplied separately in `Tesseroid-LineIntegration-data-0.1.0-rc.1.zip`. Extract that archive **into the repository root** to restore its `Works_Zhangetal2026/...` paths. Generated figures, executables, MEX binaries, manuscript drafts and reference-document copies are omitted. External ETOPO/VIM source models and map assets are listed in the reproduction guide.

## Numerical scope

The homogeneous-tesseroid identity is exact; the implementation uses finite-precision quadrature and special branches. The shell method adds interpolation and near-region source approximations. Quadrature tolerance does not bound those modelling errors. Exterior, strict-interior, exact-pole and upper-face cases follow different paths; some material-boundary values are deliberately rejected. Read [METHOD.md](docs/METHOD.md) before using high-order derivatives near a source.

The two header directories contain distinct kernel variants sharing the `Tesseroid` namespace. Use one variant per translation unit. Do not include both variants or combine their definitions in one linked program. Separate executables are provided for ordinary-TLI and TLI-MMSI validation.

## Citation and license

Use [CITATION.cff](CITATION.cff) for software authorship and the manuscript title above for the associated research. Record the final publication DOI and archive identifier when they are assigned. Software is distributed under the [MIT License](LICENSE). Third-party inputs retain their own terms; this license does not grant rights to external datasets or map assets.

## Delivery synchronization

This public repository is the canonical release record. Every later change to the source code, supplementary PDFs, figures, data documentation, or README files must be rebuilt where applicable, committed, and pushed to the `main` branch before delivery. A local-only change or a package that differs from the current GitHub contents is not a release version.
