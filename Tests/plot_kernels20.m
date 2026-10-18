%% Twenty cancellation-resistant one-dimensional tesseroid kernels
% Edit only the USER PARAMETERS section for a new prism or observation point.
% Local Cartesian axes: x = north, y = east, z = radial outward.

clc; close all;

%% USER PARAMETERS
prism.lat1_deg = 30.0;
prism.lat2_deg = 50.0;
prism.lon1_deg = 110.0;
prism.lon2_deg = 130.0;
prism.r1_km    = 6341.2;
prism.r2_km    = 6371.2;
prism.density  = 1.0;

point.lat_deg = 40.0;
point.lon_deg = 129.9;
point.r_km    = 6372.2;

nSamples = 1601;             % Number of samples on each paired line kernel
usePhysicalScale = false;   % true: multiply by G*rho*dr*dlat*dlon/4
G = 6.67430e-20;            % km^3 kg^-1 s^-2; set G=1 for normalized tests
normalizeEachPanel = false; % true: divide both curves by their common max

% Optional batch overrides used by export_near_boundary_kernels.m.
if exist('kernelOverrides','var')
    if isfield(kernelOverrides,'prism')
        f = fieldnames(kernelOverrides.prism);
        for ii=1:numel(f), prism.(f{ii})=kernelOverrides.prism.(f{ii}); end
    end
    if isfield(kernelOverrides,'point')
        f = fieldnames(kernelOverrides.point);
        for ii=1:numel(f), point.(f{ii})=kernelOverrides.point.(f{ii}); end
    end
    if isfield(kernelOverrides,'nSamples'), nSamples=kernelOverrides.nSamples; end
    if isfield(kernelOverrides,'usePhysicalScale'), usePhysicalScale=kernelOverrides.usePhysicalScale; end
    if isfield(kernelOverrides,'normalizeEachPanel'), normalizeEachPanel=kernelOverrides.normalizeEachPanel; end
    if isfield(kernelOverrides,'G'), G=kernelOverrides.G; end
end

%% Geometry and sampling
p.lat1 = deg2rad(prism.lat1_deg); p.lat2 = deg2rad(prism.lat2_deg);
p.lon1 = deg2rad(prism.lon1_deg); p.lon2 = deg2rad(prism.lon2_deg);
p.r1 = prism.r1_km; p.r2 = prism.r2_km; p.density = prism.density;
q.lat = deg2rad(point.lat_deg); q.lon = deg2rad(point.lon_deg); q.r = point.r_km;

assert(p.lat2 > p.lat1 && p.lon2 > p.lon1 && p.r2 > p.r1 && p.r1 > 0, ...
    'The prism bounds must be strictly increasing.');
assert(q.r > 0, 'The observation radius must be positive.');
assert(nSamples >= 9 && mod(nSamples,2) == 1, 'Use an odd nSamples >= 9.');

zeta = linspace(-1,1,nSamples);
lonMid = (p.lon1+p.lon2)/2; lonHalf = (p.lon2-p.lon1)/2;
latMid = (p.lat1+p.lat2)/2; latHalf = (p.lat2-p.lat1)/2;

Kphi = zeros(20,nSamples);  % latitude-boundary pair, integrated in longitude
Klam = zeros(20,nSamples);  % longitude-boundary pair, integrated in latitude

fprintf('Sampling %d points for two paired kernels...\n',nSamples);
tic;
for j = 1:nSamples
    Kphi(:,j) = pairedKernel(p,q,'latitude', lonMid+lonHalf*zeta(j));
    Klam(:,j) = pairedKernel(p,q,'longitude',latMid+latHalf*zeta(j));
end
elapsed = toc;
fprintf('Kernel sampling finished in %.3f s.\n',elapsed);

if usePhysicalScale
    commonScale = G*p.density*(p.r2-p.r1)*(p.lat2-p.lat1)*(p.lon2-p.lon1)/4;
    Kphi = commonScale*Kphi;
    Klam = commonScale*Klam;
end

%% 4 x 5 plot
names = {'V','V_x','V_y','V_z', ...
    'V_{xx}','V_{xy}','V_{xz}','V_{yy}','V_{yz}','V_{zz}', ...
    'V_{xxx}','V_{xxy}','V_{xxz}','V_{xyy}','V_{xyz}', ...
    'V_{xzz}','V_{yyy}','V_{yyz}','V_{yzz}','V_{zzz}'};

figureVisible = 'on';
if exist('kernelFigureVisible','var'), figureVisible=kernelFigureVisible; end
figure('Color','w','Name','Twenty tesseroid line kernels','Visible',figureVisible);
tl = tiledlayout(4,5,'TileSpacing','compact','Padding','compact');
for k = 1:20
    nexttile;
    y1 = Kphi(k,:); y2 = Klam(k,:);
    if normalizeEachPanel
        s = max([abs(y1),abs(y2)]);
        if s > 0, y1 = y1/s; y2 = y2/s; end
    end
    A = min([y1,y2]);
    B = max([y1,y2]);
    C = B-A;
    if C == 0
        C = max(1,abs(A));
    end
    plot(zeta,y1,'LineWidth',1.15); hold on;
    plot(zeta,y2,'LineWidth',1.15);
    yline(0,':','Color',[0.45 0.45 0.45]);
    xlim([-1 1]); ylim([A-0.1*C,B+0.1*C]); grid on; box on;
    title(['$',names{k},'$'],'Interpreter','latex');
    if k > 15, xlabel('$\xi,\eta$','Interpreter','latex'); end
    if any(k == [1 6 11 16]), ylabel('kernel'); end
    if k == 1
        legend('$K_\varphi$','$K_\lambda$','Interpreter','latex', ...
            'Location','best','FontSize',8);
    end
end
sgtitle(sprintf(['Paired one-dimensional kernels: prism [%g,%g]^\\circ x ', ...
    '[%g,%g]^\\circ, r=[%g,%g] km; point=(%g^\\circ,%g^\\circ,%g km)'], ...
    prism.lat1_deg,prism.lat2_deg,prism.lon1_deg,prism.lon2_deg, ...
    prism.r1_km,prism.r2_km,point.lat_deg,point.lon_deg,point.r_km));

if exist('kernelExportFile','var') && ~isempty(kernelExportFile)
    exportgraphics(gcf,kernelExportFile,'Resolution',240);
    fprintf('Exported: %s\n',kernelExportFile);
end

clear kernelOverrides kernelExportFile kernelFigureVisible figureVisible

%% Local functions: normalized third-order Taylor jets
function out = pairedKernel(p,point,kind,t)
    [x,r,n] = observationJet(point); %#ok<ASGLU>
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
    if ~(v>0), error('A logarithm argument is non-positive; check for a boundary singularity.'); end
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
