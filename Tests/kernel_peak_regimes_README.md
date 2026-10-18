# TLI kernel peak-regime experiment

This experiment uses the production C++ paired-boundary kernel and adaptive
Gauss--Kronrod implementation. MATLAB is retained only for inspecting kernel
shapes through `tli_paired_kernel20.m` and `plot_kernels20.m`.

## Build and run

From this directory on Windows:

```bat
build_kernel_peak_regimes.bat
kernel_peak_regimes.exe smoke results\kernel_peak_regimes\kernel_peak_regimes_cpp_smoke.csv
kernel_peak_regimes.exe position results\kernel_peak_regimes\kernel_peak_regimes_cpp_position.csv
```

Available modes are:

- `smoke`: four derivative orders, two heights, three angular gaps, and both
  paired boundary directions.
- `position`: the narrow Vzzz kernel at 19 along-line peak positions. This
  directly tests whether quadrature nodes can miss a peak.
- `quick` and `full`: broader parameter sweeps intended for unattended runs.

Each case records the shortest angular and combined spatial distance, predicted
geometric width, measured FWHM and support width, asymmetry, L1 area, reference
stability, actual integration error, embedded error estimate, node count,
interval count, convergence flag, and elapsed time.

The reference is a locally refined composite GL32 result; the GL32--GL24
difference is stored as `reference_error`. The tested methods are fixed GL6,
GL10, GL20, and GL30; unsplit adaptive GK15; and peak-localized split GK21.

## Preliminary result (smoke and position modes)

The narrow Vzzz position sweep has FWHM/interval ratios of approximately
8.5e-6 to 1.6e-5. GL6--GL30 miss all 38 tested peaks almost completely. The
ordinary adaptive GK15 did not silently accept a wrong result: it detected
unresolved structure in every position case and reached the configured 99,975
function evaluations without convergence. Thus the demonstrated failure mode
for this kernel is fixed-rule peak loss and unsplit-adaptive cost explosion.
Silent adaptive peak loss remains a hypothesis to test for other kernels and
geometries, not a conclusion from this data set.

Simple interval splitting alone is also not sufficient for the most extreme
third-derivative cases because roundoff in the near-singular kernel can drive
the embedded estimator to its evaluation ceiling. The production estimator's
cosine/tanh--sinh transforms should therefore be compared as an additional
method before setting final switching thresholds.
