function experiment = test1_geometry(heightAboveTopMetres)
%TEST1_GEOMETRY Non-symmetric source and an angularly off-centre point.

arguments
    heightAboveTopMetres (1,1) double {mustBePositive}
end

prism.lat1 = 8.3;
prism.lat2 = 11.1;
prism.lon1 = 19.7;
prism.lon2 = 23.4;
prism.r1 = 6.3412e6;
prism.r2 = 6.3712e6;
prism.density = 1;

angularCentre = [0.5*(prism.lat1+prism.lat2), ...
    0.5*(prism.lon1+prism.lon2)];
point = [angularCentre+[1,1],prism.r2+heightAboveTopMetres];
assert(point(1)>prism.lat1 && point(1)<prism.lat2);
assert(point(2)>prism.lon1 && point(2)<prism.lon2);

experiment.prism = prism;
experiment.angularCentreDegrees = angularCentre;
experiment.angularOffsetDegrees = [1,1];
experiment.point = point;
experiment.heightAboveTopMetres = heightAboveTopMetres;
experiment.axes = "x=north, y=east, z=radially outward at the observation point";
experiment.componentOrder = ...
    "V,Vx,Vy,Vz,Vxx,Vxy,Vxz,Vyy,Vyz,Vzz," + ...
    "Vxxx,Vxxy,Vxxz,Vxyy,Vxyz,Vxzz,Vyyy,Vyyz,Vyzz,Vzzz";
end
