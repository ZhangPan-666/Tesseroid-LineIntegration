function value = t1s_vpa2d_component(prism,point,componentIndex, ...
    workingDigits,quadratureDigits)
%T1S_VPA2D_COMPONENT Analytically radial-integrated 2-D VPA integral.

oldDigits = digits(workingDigits);
restoreDigits = onCleanup(@()digits(oldDigits));
opts = t1s_vpa_options(workingDigits,quadratureDigits);
g = t1s_normalized_geometry(prism,point,workingDigits);
phi = sym("phi","real");
lambda = sym("lambda","real");
kernel = t1s_kernel2d_component(g,phi,lambda,componentIndex);

latitudinal = vpaintegral(kernel,phi,[g.lat1,g.lat2], ...
    "RelTol",opts.relTol,"AbsTol",opts.absTol);
value = vpaintegral(latitudinal,lambda,[g.lon1,g.lon2], ...
    "RelTol",opts.relTol,"AbsTol",opts.absTol);
value = vpa(value,workingDigits);
end
