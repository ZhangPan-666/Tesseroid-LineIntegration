# Method and implementation scope

This description follows the supplied Zhang–Du manuscript and the packaged source. The manuscript and supplementary derivations are not redistributed in this software archive.

## Homogeneous tesseroid TLI

For positive Newtonian potential, start with the divergence field `F0=(y-x)/(2|y-x|)`, whose source divergence is `1/|y-x|`. A curl correction leaves that divergence unchanged and makes normal flux constant on a source sphere. The spherical-face contribution becomes explicit. Analytical integration in source radius removes the radial variable from the four lateral faces, leaving two oriented pairs of boundary-curve integrals.

In the manuscript notation,

`V = (G rho / 2) [C0 + I_latitude + I_longitude]`,

where `C0=(r2² m(r2)-r1² m(r1)) Δlongitude (sin(lat2)-sin(lat1))`. The sphere-mean coefficient is `m(a)=1-r²/(3a²)` for `r<a`, `2/3` at equality and `2a/(3r)` for `r>a`.

The production kernel uses third-order multivariate Taylor jets for Cartesian derivatives. Rationalized radial differences, endpoint branches, `log1p` and local source-parameter expansions reduce cancellation. Opposite boundaries are paired before accumulation, with compensated summation. The local basis is fixed while differentiating.

Quadrature splits at an analytical closest-angular-projection candidate on the relevant boundary line. Cosine/sine maps cluster nodes near that projection. Adaptive GK21 is selected for narrow nearby kernels; GK15 is used otherwise. Each transformed piece is capped at 10,000 evaluations and 500 intervals, subject to any lower caller cap. Relative tolerance has a binary64 floor of `500*epsilon`; a kernel-magnitude roundoff floor also enters acceptance. Requesting `1e-16` does not establish 16-digit accuracy.

## Boundary and pole handling

- Strict exterior and strict interior points use the combined reduced representation, including the Poisson contribution in matter.
- On the open upper radial face, the implementation returns an exterior one-sided Richardson-extrapolated value. This is a separate numerical operation.
- On the lower radial material face and latitude/longitude faces, edges or corners, ordinary evaluation rejects an undefined boundary convention.
- At exact poles, if the observation is radially or angularly separated from the source, a special Cartesian calculation analytically integrates longitude and applies a double-exponential product rule in the other two variables. This branch precedes material-boundary classification. Other pole cases follow the remaining ordinary logic.

These distinctions matter for high derivatives. A finite quadrature value at a material discontinuity is not sufficient to define the intended one-sided tensor. Symmetry/trace identities also do not independently validate every component.

## Gridded-shell TLI-MMSI

The shell algorithm reconstructs upper/lower radii and material properties on a regular spherical grid. It uses three numerical regions around each observation direction:

1. `H0`, one observation-centred source-grid cell, receives high-accuracy TLI.
2. `H1\H0`, extending at least one degree on each side, aggregates shared oriented boundaries, then integrates interpolation panels with GL2/GL4 and clustered GL10 near the projection.
3. Four far rectangles outside `H1` use adaptive two-dimensional radially reduced surface kernels. The embedded rule operates on 14×14 node arrays, with an evaluation-budget parameter of 500 per rectangle and component.

For cells wider than two degrees, the near region reduces to the central TLI cell plus the far calculation. Component masks avoid unnecessary kernels for Vz-only gravity and magnetic-vector-only calculations. The far evaluator internally obtains `Integral2Info` but the public shell interface does not propagate those convergence/budget warnings; requested tolerances are therefore not an end-to-end shell error guarantee.

The central-cell wrapper floors both absolute and relative requests at `1e-10`, in addition to the standalone floating-point controls. Effective settings thus differ between numerical regions.

Shell interpolation and the homogeneous near-cell replacement introduce source-model approximations. The shell output is not the exact response of an arbitrary collection of original independent tesseroids. Compare methods using an explicitly shared source representation before making pointwise accuracy or speed claims.

## Distinct header variants

The `C++Codes` and `TLI_MMSI` production headers were retained byte-for-byte. The shell kernel includes a short thin-layer radial-series branch absent from ordinary TLI. Identical-looking names do not imply identical implementation. Keep variants in separate executables; the release tests do so. Duplicate ordinary headers retained in the standalone cluster/runtime folders have matching original hashes, documented by `SOURCE_INVENTORY.csv`.
