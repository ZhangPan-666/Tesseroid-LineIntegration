# Single-tesseroid MPI runtime benchmark

Build `Test_RuntimeTest` with root CMake, MPI and research ON, or run `make CXX=mpicxx` here. Run under MPI with this folder as working directory. `run.pbs` is a generic SLURM example; adjust resource requests for the cluster.

The benchmark contains 160,801 observations and per-component warmups/repetitions. Review its constants before allocating a full run. Each MPI process is single-threaded; no OpenMP flag is added.

`Analyze_RuntimeTest1.m` reads `RuntimeTest_LaplacePoisson_HeightPlus1km.csv` (in the data attachment), summarizes per-component/derivative-rank timings and projection regions, and writes `results/`. `Summarize_RuntimeTest.m` is a compatibility entry point. `Figures/Figure_RuntimeTest_1km.m` plots archived timings. Measurements depend on compiler, hardware and settings and are not release-check timings.

See [reproduction](../../docs/REPRODUCIBILITY.md).
