# Build and run

## C++

Use CMake 3.18+ and a C++17 compiler. Configure at the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
```

For Visual Studio, use a Developer PowerShell or let CMake select the installed generator. `--config Release` selects optimized code on multi-configuration generators. With GCC/Clang it is selected by `CMAKE_BUILD_TYPE`.

Optional switches:

| Switch | Default | Adds |
| --- | --- | --- |
| `BUILD_TESTING` | ON | Independent volume reference, shell reference and smoke checks |
| `TLI_BUILD_OPENMP` | OFF | Ordinary TLI OpenMP drivers; requires OpenMP |
| `TLI_BUILD_MPI` | OFF | Ordinary TLI and TLI-MMSI MPI drivers; requires MPI C++ development files |
| `TLI_BUILD_RESEARCH` | OFF | Long research programs and, with MPI, shell/runtime benchmark programs |

```sh
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DTLI_BUILD_MPI=ON -DTLI_BUILD_RESEARCH=ON
cmake --build build-mpi --config Release --parallel 2
```

Every MPI executable links MPI only. OpenMP targets are separate even when both build switches are enabled. Do not add global OpenMP flags to MPI builds. Example cluster use, from the shell-validation working directory:

```sh
mpirun -np 4 /absolute/path/to/build-mpi/shell_validation 60
```

The shell argument is an observation sampling stride; inspect its printed usage before a full run. The `run.pbs` files contain SLURM directives despite their historical suffix. Edit cluster resource requests to suit your system; they default to a small task count in this package.

The legacy Makefiles also support `make CXX=mpicxx`. For an Intel MPI installation, pass its C++ compiler wrapper through `CXX`; no Intel-specific OpenMP flags are inserted by these Makefiles.

After building both OpenMP and MPI targets, the optional standard-library Python check compares the batch entry points (use `build/Release` with Visual Studio):

```sh
python tools/check_driver_modes.py --serial-dir build --mpi-dir build --scratch driver-checks
```

It launches two local ranks through `mpiexec`, writes small synthetic parameter files in the dedicated scratch folder, and compares serial/OpenMP/MPI and full/component-selective outputs. Add your MPI runtime to the executable search path first if required.

## Work file drivers

Ordinary drivers read `TFG.ForPar` or `TFM.ForPar` and write `Result_TFG.dat` or `Result_TFM.dat` in the **current directory**. Run the matching `C++Codes/Work1_Set*Parameters.m` generator there first. Shell drivers read `MSIG.ForPar` or `MSIM.ForPar` and write `ResultG.dat` or `ResultM.dat`. Use distinct working directories for different runs to avoid overwriting results.

These retained batch programs assume correctly formed input files. The direct C++ API and MEX gateway are preferable for new integrations. Binary files use native C++ `int32`/binary64 records; archived files were generated on little-endian systems. See [API.md](API.md).

## MATLAB MEX

Add `Functions/` to the MATLAB path. Configure a supported C++ compiler using `mex -setup C++`, then run:

```matlab
build_Tesseroid_LineIntegration_mex       % OpenMP by default on Windows/Linux
build_Tesseroid_LineIntegration_mex(false) % serial MEX; also default on macOS
run('examples/demo_matlab.m')
```

MATLAB R2025b is the research environment. Other releases have not been tested. OpenMP compiler support is required only for the `P` wrappers. On macOS the default serial build uses Apple Clang; a custom OpenMP toolchain needs appropriate compiler and linker configuration.

An optional Intel/Windows build script remains available. It requires `TLI_MATLAB_ROOT` to name the MATLAB installation and optionally `TLI_ONEAPI_ROOT` (default: `C:\Program Files (x86)\Intel\oneAPI`). Run `Functions/build_mex_intel.bat` after setting those values. Its output is a local build artifact and is excluded from the source ZIP.

The MEX gateway requires all 20/10 internal outputs; the MATLAB wrappers already provide them even if the caller requests fewer. Use the wrappers rather than calling the gateway directly.

## Verification status

Local compiler versions, test results, and untested platforms are recorded in [VALIDATION.md](VALIDATION.md). GitHub Actions definitions are provided; they become executed evidence only after the repository is pushed and the jobs pass.
