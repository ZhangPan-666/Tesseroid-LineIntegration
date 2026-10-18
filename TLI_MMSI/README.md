# TLI-MMSI hybrid method

The gravity and magnetization programs use one consistent source model with
three numerical regions.

- H0 is the observation-centred tesseroid with the source-grid dimensions. It
  is evaluated by high-accuracy adaptive TLI.
- H1 is the observation-centred rectangle extending at least one degree in
  longitude and latitude. Its part outside H0 is evaluated by aggregated TLI
  on source-interpolation panels (GL2/GL4, with cosine-clustered GL10 beside
  the projection).
- The four rectangles outside H1 are evaluated by adaptive two-dimensional
  surface integration without peak finding. Each rectangle is limited to 500
  evaluations of the 14 by 14 embedded quadrature rule.

If a source cell is wider than two degrees, H0 already contains the nominal
one-degree near field and the method reduces to H0 TLI plus two-dimensional
far-field integration. Component masks are respected, so single-component and
magnetic calculations evaluate only the surface kernels they require.

The two-dimensional integration dependencies are `_integral2.h`,
`_sphinterp.h`, and `_tesseroid_surfacekernel.h`. Fixed 14 by 14 quadrature
arrays use the C++ standard library, so the programs have no Eigen dependency.

## Work2 component-selective calculations

`Work2_MassShellIntegrationVzMPI` evaluates only the canonical `Vz` component.
`Work2_MagnetizationShellIntegrationVectorMPI` evaluates only the nine Hessian
components required to assemble magnetic `Vx`, `Vy`, and `Vz`; it does not
evaluate the magnetic potential or magnetic-gradient outputs. The magnetic
Work2 implementation also reuses its transformed-magnetization work arrays and
precomputed source-angle trigonometric tables for all observations assigned to
an MPI rank.

Work1 and Work2 timers both include source/observation broadcasts, numerical
integration, and result gathers. They synchronize before and after the timed
region and report the maximum elapsed time over all MPI ranks. Result-file I/O
is outside the printed integration time.
