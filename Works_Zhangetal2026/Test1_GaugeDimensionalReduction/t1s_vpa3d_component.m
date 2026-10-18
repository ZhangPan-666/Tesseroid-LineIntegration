function value = t1s_vpa3d_component(prism,point,componentIndex, ...
    workingDigits,quadratureDigits)
%T1S_VPA3D_COMPONENT Original 3-D VPA derivative-kernel integral.

oldDigits = digits(workingDigits);
restoreDigits = onCleanup(@()digits(oldDigits));
opts = t1s_vpa_options(workingDigits,quadratureDigits);
g = t1s_normalized_geometry(prism,point,workingDigits);

a = sym("a","real");
phi = sym("phi","real");
lambda = sym("lambda","real");
[tx,ty,cosine] = t1s_geometry_terms(phi,lambda,g);
dx = a*tx;
dy = a*ty;
dz = a*cosine-g.r;
distance = sqrt(dx^2+dy^2+dz^2);
kernels = t1s_kernel20(dx,dy,dz,distance,a^2*cos(phi));
kernel = kernels(componentIndex);

radial = vpaintegral(kernel,a,[g.r1,g.r2], ...
    "RelTol",opts.relTol,"AbsTol",opts.absTol);
latitudinal = vpaintegral(radial,phi,[g.lat1,g.lat2], ...
    "RelTol",opts.relTol,"AbsTol",opts.absTol);
value = vpaintegral(latitudinal,lambda,[g.lon1,g.lon2], ...
    "RelTol",opts.relTol,"AbsTol",opts.absTol);
value = vpa(value,workingDigits);
end

