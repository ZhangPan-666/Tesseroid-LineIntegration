# 0.1.0-rc.1 — publication candidate

This candidate packages homogeneous-tesseroid line integration for potential and Cartesian derivatives through third order, MATLAB serial/OpenMP interfaces, and hybrid TLI-MMSI mass/magnetization shell drivers.

It adds a CMake build, separate pure-MPI and OpenMP targets, small examples, independent numerical checks, bilingual documentation, API/unit guidance, MIT licensing and reproducibility manifests. Production numerical headers are preserved from the supplied snapshot. The MEX gateway now preserves validation error identifiers.

Local Windows checks passed: the three CTest tests, ordinary serial/OpenMP/two-rank MPI comparisons, shell component-selection comparisons, MATLAB gravity/magnetic interface checks, and archived validation/runtime analysis. Linux/macOS CI is configured but has not yet run. See `docs/VALIDATION.md` for scope and observed errors.

Attach both source and data ZIPs plus `SHA256SUMS.txt`. The data attachment restores archived outputs when extracted into the source root. It excludes external global source models and map assets. Historical full-grid arrays are not relabelled as fresh calculations of this candidate; exact source-version correspondence and some external reproduction dependencies remain unresolved.

The associated manuscript is by Pan Zhang and Jinsong Du, “An exact line-integral formulation for high-accuracy evaluation of tesseroid potentials and their derivatives through third order.” Add the final manuscript/software DOI and repository URL when assigned.
