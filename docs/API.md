# API, coordinates and units

## Ordinary C++ API

Include `_TLI_estimate.h` from `C++Codes/`.

```cpp
Tesseroid::Prism p{lat1,lat2,lon1,lon2,r1,r2,density};
Tesseroid::Point x{latitude,longitude,radius};
Tesseroid::EstimateOptions options;
auto all = Tesseroid::estimate_all(p,x,options);
Tesseroid::EstimateInfo info;
double vz = Tesseroid::estimate(Tesseroid::Component::Vz,p,x,options,&info);
```

Angles are **geocentric radians**. Require `0 < r1 < r2`, `radius > 0`, `-pi/2 <= lat1 < lat2 <= pi/2`, `-pi <= lon1 < lon2 <= pi`, and `-pi <= observation longitude < pi`. Split antimeridian-crossing source cells first. The MATLAB and Work adapters normalize observation longitude; the direct API validates it. Density may be signed; geometry and inputs must be finite.

Outputs use a fixed local Cartesian frame at the observation: x north, y east, z radially outward. The potential is positive `G integral(rho/d) dvolume`, with gravity `grad(V)`. Thus positive mass normally produces negative exterior `Vz`; the material-interior trace is `Vxx + Vyy + Vzz = -4*pi*G*rho`.

At an exact pole the implementation fixes the horizontal axes using the zero-longitude reference meridian, regardless of the input pole longitude.

`EstimateInfo` provides `value`, estimated absolute quadrature `error`, `absolute` kernel magnitude, `integral_evaluations`, `intervals`, `converged` and `milliseconds`. The estimate is an adaptive stopping indicator, not a rigorous total-error bound. Always inspect `converged` for requested precision. The Work adapters discard these diagnostics; magnetic MEX checks finite values but does not expose the underlying convergence flags.

## Component mapping

Positions below are one-based. Alias names are equal mixed partial derivatives in the fixed Cartesian frame; the **permutation also matters**.

| Position | Canonical `estimate_all` | MATLAB / Work / public shell output |
| --- | --- | --- |
| 1 | V | V |
| 2 | Vx | Vx |
| 3 | Vy | Vy |
| 4 | Vz | Vz |
| 5 | Vxx | Vxx |
| 6 | Vxy | Vxy |
| 7 | Vxz | Vyy |
| 8 | Vyy | Vzx (= Vxz) |
| 9 | Vyz | Vzy (= Vyz) |
| 10 | Vzz | Vzz |
| 11 | Vxxx | Vxxx |
| 12 | Vxxy | Vxxy |
| 13 | Vxxz | Vxxz |
| 14 | Vxyy | Vxyz |
| 15 | Vxyz | Vyyx (= Vxyy) |
| 16 | Vxzz | Vyyy |
| 17 | Vyyy | Vyyz |
| 18 | Vyyz | Vzzx (= Vxzz) |
| 19 | Vyzz | Vzzy (= Vyzz) |
| 20 | Vzzz | Vzzz |

The ten magnetic outputs are `V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz`. They are magnetization-weighted gravity derivative combinations. Magnetization vectors are supplied in source-centred north/east/up coordinates and rotated into the observation frame.

## Units

For the direct API, the default `gravitational_constant` is `6.67430e-11` in SI. Metres and kg/m³ give:

| Derivative order | SI unit |
| --- | --- |
| Potential | m²/s² |
| First | m/s² |
| Second | 1/s² |
| Third | 1/(m s²) |

Compatibility interfaces (`Tesseroid_GravityEstimate*`, `Tesseroid_MagneticEstimate*`, MATLAB, and shell adapters) use **G=1** and degrees. They neither insert SI G nor automatically convert km. For raw gravity values with density numerically in kg/m³ and coordinates in km, order `n` converts to SI by multiplying by `G_SI * 1000^(2-n)`. This follows from the volume and derivative length scaling. Multiply SI acceleration by `1e5` for mGal, Hessian by `1e9` for Eötvös, and third derivatives by `1e12` for mE/km.

The magnetic adapters do not include `mu0/(4*pi)`. With coordinates in metres and magnetization in A/m, multiply raw potential, field and field-gradient quantities by the chosen SI magnetic constant `mu0/(4*pi)`. With km coordinates, first convert their respective raw scales by `1000`, `1` and `0.001`. For the legacy nominal `mu0/(4*pi)=1e-7`, these give potential ×0.1 (microtesla metre), field ×100 (nT), and gradient ×100 (nT/km), matching the retained plotting conventions. Sign and output definitions follow the implemented derivative combinations, so document your chosen magnetic-potential convention when exporting results.

## MATLAB matrices

Gravity model: N×7 `[lat1 lat2 lon1 lon2 r1 r2 density]`.
Magnetic model: N×9 `[lat1 lat2 lon1 lon2 r1 r2 Mx My Mz]`.
Observation matrix: M×3 `[latitude longitude radius]`.
Matrices are full, real doubles. `AT` and `RT` must be finite nonnegative scalars. Angles are degrees. All radii must use the same length unit.

## Shell and binary layouts

Shell fields are latitude-fastest: MATLAB matrices are `NumLat × NumLon`, and C++ indexes `i_lon * NumLat + i_lat`. Grid limits describe **cell centres**. Use the shell API for its intended global regular grid; a small regional grid smoke check does not establish regional-shell modelling accuracy.

Ordinary parameter files contain doubles `AT,RT`, int32 model count, seven (gravity) or nine (magnetic) model-column arrays, int32 observation count, then `longitude,latitude,radius` arrays. Ordinary result files contain int32 observation count, those three coordinate arrays, then all output-component arrays in the public order above.

Shell parameter files begin with doubles `AT,RT,lon_first,lon_last,lat_first,lat_last,dlon,dlat`, int32 `NumLon,NumLat`, then `R2,R1,rho` (gravity) or `R2,R1,Mx,My,Mz` (magnetic), int32 observation count, and the three observation arrays. Retained shell result files use the ordinary result layout with 20 gravity or 10 magnetic arrays; elapsed time is printed separately. Component-selective Work2 files keep the Work1 output shape: inspect the driver for unused-field zero/NaN fill values.
