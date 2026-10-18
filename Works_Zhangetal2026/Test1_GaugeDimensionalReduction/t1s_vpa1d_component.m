function value = t1s_vpa1d_component(prism,point,componentIndex, ...
    workingDigits,quadratureDigits)
%T1S_VPA1D_COMPONENT Differentiate the gauge-reduced scalar line formula.

oldDigits = digits(workingDigits);
restoreDigits = onCleanup(@()digits(oldDigits));
opts = t1s_vpa_options(workingDigits,quadratureDigits);
g = t1s_normalized_geometry(prism,point,workingDigits);
meta = t1s_components();
exponent = meta.exponents(componentIndex,:);

X=sym("X","real"); Y=sym("Y","real"); Z=sym("Z","real");
t=sym("t","real");
n0 = [cos(g.lat)*cos(g.lon);cos(g.lat)*sin(g.lon);sin(g.lat)];
ex = [-sin(g.lat)*cos(g.lon);-sin(g.lat)*sin(g.lon);cos(g.lat)];
ey = [-sin(g.lon);cos(g.lon);sym(0)];
x = g.r*n0+X*ex+Y*ey+Z*n0;
ro = sqrt(sum(x.^2));

branch1 = local_branch(g.r1,g.r);
branch2 = local_branch(g.r2,g.r);
m1 = t1s_gauge_mean(g.r1,ro,branch1);
m2 = t1s_gauge_mean(g.r2,ro,branch2);
omega = (g.lon2-g.lon1)*(sin(g.lat2)-sin(g.lat1));
c0 = (g.r2^2*m2-g.r1^2*m1)*omega;
c0d = local_derivative(c0,X,Y,Z,exponent);
c0d = subs(c0d,[X,Y,Z],[0,0,0]);

epsilon = [-1,1];
latBounds = [g.lat1,g.lat2];
lonBounds = [g.lon1,g.lon2];
edgeSum = sym(0);

for k=1:2
    phi = latBounds(k);
    u = [cos(phi)*cos(t);cos(phi)*sin(t);sin(phi)];
    ephi = [-sin(phi)*cos(t);-sin(phi)*sin(t);cos(phi)];
    c = dot(x,u)/ro;
    cphi = dot(x,ephi)/ro;
    w = local_w(g.r1,g.r2,ro,c,branch1,branch2);
    integrand = epsilon(k)*cos(phi)*cphi*w;
    differentiated = local_derivative(integrand,X,Y,Z,exponent);
    differentiated = subs(differentiated,[X,Y,Z],[0,0,0]);
    edgeSum = edgeSum+vpaintegral(differentiated,t, ...
        [g.lon1,g.lon2],"RelTol",opts.relTol,"AbsTol",opts.absTol);
end

for k=1:2
    lambda = lonBounds(k);
    u = [cos(t)*cos(lambda);cos(t)*sin(lambda);sin(t)];
    elambda = [-sin(lambda);cos(lambda);sym(0)];
    c = dot(x,u)/ro;
    clambda = dot(x,elambda)/ro;
    w = local_w(g.r1,g.r2,ro,c,branch1,branch2);
    integrand = epsilon(k)*clambda*w;
    differentiated = local_derivative(integrand,X,Y,Z,exponent);
    differentiated = subs(differentiated,[X,Y,Z],[0,0,0]);
    edgeSum = edgeSum+vpaintegral(differentiated,t, ...
        [g.lat1,g.lat2],"RelTol",opts.relTol,"AbsTol",opts.absTol);
end

value = vpa((c0d+edgeSum)/2,workingDigits);
end

function branch = local_branch(a,r)
if isAlways(a<r)
    branch="source_inside_observer";
elseif isAlways(a>r)
    branch="source_outside_observer";
else
    error("Observation radius equals a source radial endpoint.");
end
end

function w = local_w(r1,r2,r,c,branch1,branch2)
w = t1s_gauge_q(r2,r,c,branch2)- ...
    t1s_gauge_q(r1,r,c,branch1)-t1s_radial_h(r1,r2,r,c);
end

function value = local_derivative(value,X,Y,Z,exponent)
for k=1:exponent(1), value=diff(value,X); end
for k=1:exponent(2), value=diff(value,Y); end
for k=1:exponent(3), value=diff(value,Z); end
end

