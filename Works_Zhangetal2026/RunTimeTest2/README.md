# Global topographic mass-shell benchmark

Archived result arrays are in the separate data attachment. `RunTimeRecord.m` retains historical measured times; these were not remeasured during release packaging.

To regenerate inputs, configure `TLI_ETOPO_DIR` with the resampled ETOPO tables and run the `RunTimeTest_Gravity.m` generator. Run the matching pure-MPI shell executable from a working directory containing the generated `MSIG.ForPar` file. Keep separate output directories for all-components and component-selective runs.

The original map plotting script additionally requires `TLI_GEOCONTOUR_DIR`, MATLAB Mapping Toolbox and the exact helper/shapefile/MAT assets it names. These third-party inputs were absent from the supplied directory. No data-source download or redistribution right is inferred.

Current TLI-MMSI uses observation-centred near reconstruction and 2D far-field integration. Historical notes about a pure grid-line model or exploratory speedup files missing from the archive have been removed. Old arrays are not claimed to reproduce the current source snapshot without a full rerun. See [reproduction](../../docs/REPRODUCIBILITY.md) and [units](../../docs/API.md).
