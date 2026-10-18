# Unified validation analysis

`RunValidationAnalysis.m` is the single statistical audit for manuscript
Tests No. 1--11. It reads the archived Test1--Test4 result files without
modifying them and reports, for every relevant error sample:

- sample count and nonfinite count;
- minimum, maximum, arithmetic mean, and RMSE;
- empirical P50, P90, and P99.

Run with MATLAB R2025b:

```matlab
run('RunValidationAnalysis.m')
```

Generated files are written to `Analysis/results/`:

- `validation_audit.log`: human-readable console audit, including every
  manuscript-level maximum residual and timing claim;
- `validation_component_statistics.csv`: full component/identity statistics;
- `validation_experiment_statistics.csv`: pooled experiment summaries;
- `manuscript_claims.csv`: compact values used by the main text;
- `validation_statistics_tables.tex`: reproducible SI Tables S1--S11.

The script treats relative errors near analytical zero crossings as valid
data rather than silently removing them. The percentile columns make their
local character visible, while the direct-value figures and absolute-error
tables preserve the physical interpretation.
