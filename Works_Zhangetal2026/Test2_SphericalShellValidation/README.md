# Complete spherical-shell validation

The driver defines a 1-degree global grid, radii 6341.2 and 6371.2 km, G=1 and density=1. It compares 20 components with closed forms outside the shell (+1 km), in the cavity (-31 km) and in matter (-15 km).

Build `shell_validation` from root CMake with `TLI_BUILD_MPI=ON` and `TLI_BUILD_RESEARCH=ON`. From this directory run `mpirun -np 4 /path/to/shell_validation 60` for a sampled calculation; stride 1 is the complete 64,800-observation grid. An optional second argument selects legacy test number 5, 6 or 7. These legacy numbers differ from the latest manuscript numbering (6–8).

Outputs overwrite `results/shell_summary.csv` and the matching `results/ResultsPan_10x10Split_<height>.dat` files. `Figures/plot_plus1km.m`, `plot_minus31km.m` and `plot_minus15km.m` reproduce maps from archived arrays supplied in the data attachment.

The old README described far-field GL6 line integration. Current source uses the hybrid 2D far field. The archived full-grid arrays have no precise source-version link and are not claimed to be new runs of the current candidate. Consult [reproduction details](../../docs/REPRODUCIBILITY.md) before citing their accuracy.
