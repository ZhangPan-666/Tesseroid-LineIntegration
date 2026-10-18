function values = tli_paired_kernel20(prism, point, kind, parameter)
%TLI_PAIRED_KERNEL20 Evaluate the 20 paired TLI boundary kernels.
%   VALUES = TLI_PAIRED_KERNEL20(PRISM,POINT,KIND,PARAMETER) returns a
%   20-by-numel(PARAMETER) matrix. KIND is "latitude" for the paired
%   latitude boundaries integrated in longitude, or "longitude" for the
%   paired longitude boundaries integrated in latitude. Angles and
%   PARAMETER are in radians; radii use one consistent length unit.
%
%   PRISM fields: lat1, lat2, lon1, lon2, r1, r2, density.
%   POINT fields: lat, lon, r.
%
%   Component order:
%   V,Vx,Vy,Vz,Vxx,Vxy,Vxz,Vyy,Vyz,Vzz,
%   Vxxx,Vxxy,Vxxz,Vxyy,Vxyz,Vxzz,Vyyy,Vyyz,Vyzz,Vzzz.

arguments
    prism struct
    point struct
    kind {mustBeMember(kind,{'latitude','longitude'})}
    parameter double
end

t = parameter(:).';
values = zeros(20,numel(t));
for i = 1:numel(t)
    values(:,i) = pairedKernel(prism,point,char(kind),t(i));
end
end

function out = pairedKernel(p,point,kind,t)
    [~,r,n] = observationJet(point);
    if strcmp(kind,'latitude')
        width = p.lat2-p.lat1;
        [um,tm] = sourceGeometry('latitude',p.lat1,t);
        [up,tp] = sourceGeometry('latitude',p.lat2,t);
        am = jscale(jdot(n,tm),cos(p.lat1));
        ap = jscale(jdot(n,tp),cos(p.lat2));
    else
        width = p.lon2-p.lon1;
        [um,tm] = sourceGeometry('longitude',p.lon1,t);
        [up,tp] = sourceGeometry('longitude',p.lon2,t);
        am = jdot(n,tm);
        ap = jdot(n,tp);
    end
    cm = jdot(n,um); cp = jdot(n,up);
    wm = wHat(p,r,cm); wp = wHat(p,r,cp);
    averageA = jscale(jadd(ap,am),0.5);
    averageW = jscale(jadd(wp,wm),0.5);
    ddW = jscale(jsub(wp,wm),1/width);
    ddA = jscale(jsub(ap,am),1/width);
    jet = jadd(jmul(averageA,ddW),jmul(averageW,ddA));
    out = jetDerivatives(jet);
end

function [x,r,n] = observationJet(p)
    er = radial(p.lat,p.lon);
    ex = north(p.lat,p.lon);
    ey = east(p.lon);
    X = jvar(1); Y = jvar(2); Z = jvar(3);
    x = cell(3,1);
    for k = 1:3
        x{k} = jadd(jconst(p.r*er(k)), ...
            jadd(jscale(X,ex(k)),jadd(jscale(Y,ey(k)),jscale(Z,er(k)))));
    end
    r = jsqrt(jadd(jmul(x{1},x{1}),jadd(jmul(x{2},x{2}),jmul(x{3},x{3}))));
    n = {jdivide(x{1},r);jdivide(x{2},r);jdivide(x{3},r)};
end

function w = wHat(p,r,c)
    dr = p.r2-p.r1;
    w = jscale(jsub(uEndpoint(p.r2,r,c),uEndpoint(p.r1,r,c)),1/dr);
end

function u = uEndpoint(a,r,c)
    u = jsub(qEndpoint(a,r,c),sourcePrimitive(a,r,c));
end

function s = sourcePrimitive(a,r,c)
    d = distanceJet(a,r,c);
    arg = jadd(jsub(jconst(a),jmul(r,c)),d);
    s = jadd(jmul(r,d),jmul(jmul(jmul(r,r),c),jlog(arg)));
end

function q = qEndpoint(a,r,c)
    d = distanceJet(a,r,c);
    r2 = jmul(r,r);
    if r(1) < a
        num = jmul(r,jadd(jsub(r2,jconst(a*a)),jscale(d,a)));
        den = jsub(jmul(jadd(jconst(a),d),jadd(jconst(a),d)),r2);
        q = jscale(jdivide(num,den),2*a/3);
    else
        num = jadd(jmul(r,d),jscale(jsub(r2,jconst(a*a)),2));
        den = jsub(jmul(jadd(r,d),jadd(r,d)),jconst(a*a));
        q = jscale(jdivide(num,den),2*a*a/3);
    end
end

function d = distanceJet(a,r,c)
    d = jsqrt(jsub(jadd(jconst(a*a),jmul(r,r)),jscale(jmul(r,c),2*a)));
end

function [u,tangent] = sourceGeometry(kind,fixed,parameter)
    if strcmp(kind,'latitude')
        u = radial(fixed,parameter);
        tangent = north(fixed,parameter);
    else
        u = radial(parameter,fixed);
        tangent = east(fixed);
    end
end

function v = radial(lat,lon)
    v = [cos(lat)*cos(lon);cos(lat)*sin(lon);sin(lat)];
end
function v = north(lat,lon)
    v = [-sin(lat)*cos(lon);-sin(lat)*sin(lon);cos(lat)];
end
function v = east(lon)
    v = [-sin(lon);cos(lon);0];
end
function z = jdot(a,b)
    z = jconst(0);
    for k = 1:3, z = jadd(z,jscale(a{k},b(k))); end
end

function d = jetDerivatives(x)
    mi = multiIndices();
    d = x.*prod(factorial(mi),2);
end
function x = jconst(c)
    x = zeros(20,1); x(1) = c;
end
function x = jvar(axis)
    x = jconst(0); x(axis+1) = 1;
end
function z = jadd(x,y), z = x+y; end
function z = jsub(x,y), z = x-y; end
function z = jscale(x,a), z = a*x; end

function z = jmul(x,y)
    persistent I J K
    if isempty(I)
        mi = multiIndices(); lookup = zeros(4,4,4);
        for k = 1:20
            lookup(mi(k,1)+1,mi(k,2)+1,mi(k,3)+1) = k;
        end
        I=[]; J=[]; K=[];
        for i=1:20
            for j=1:20
                s=mi(i,:)+mi(j,:);
                if sum(s)<=3
                    I(end+1,1)=i; J(end+1,1)=j; %#ok<AGROW>
                    K(end+1,1)=lookup(s(1)+1,s(2)+1,s(3)+1); %#ok<AGROW>
                end
            end
        end
    end
    z = accumarray(K,x(I).*y(J),[20 1],@sum,0);
end

function z = jdivide(x,y), z = jmul(x,jinverse(y)); end
function z = jinverse(x)
    v=x(1); z=jcompose(x,1/v,-1/v^2,1/v^3,-1/v^4);
end
function z = jsqrt(x)
    v=x(1); s=sqrt(v);
    z=jcompose(x,s,1/(2*s),-1/(8*v*s),1/(16*v*v*s));
end
function z = jlog(x)
    v=x(1);
    if ~(v>0)
        error('A logarithm argument is non-positive; check for a boundary singularity.');
    end
    z=jcompose(x,log(v),1/v,-1/(2*v^2),1/(3*v^3));
end
function z = jcompose(x,f0,f1,f2,f3)
    h=x; h(1)=0; h2=jmul(h,h);
    z=jadd(jconst(f0),jadd(jscale(h,f1), ...
        jadd(jscale(h2,f2),jscale(jmul(h2,h),f3))));
end
function mi = multiIndices()
    mi=[0 0 0;1 0 0;0 1 0;0 0 1;2 0 0;1 1 0;1 0 1;0 2 0;0 1 1;0 0 2; ...
        3 0 0;2 1 0;2 0 1;1 2 0;1 1 1;1 0 2;0 3 0;0 2 1;0 1 2;0 0 3];
end
