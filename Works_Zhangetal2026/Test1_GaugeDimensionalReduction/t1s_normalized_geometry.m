function g = t1s_normalized_geometry(prism,point,workingDigits)
%T1S_NORMALIZED_GEOMETRY Convert double-defined input to VPA/radians.

g.scale = vpa(max(point(3),prism.r2),workingDigits);
g.r = vpa(point(3),workingDigits)/g.scale;
g.r1 = vpa(prism.r1,workingDigits)/g.scale;
g.r2 = vpa(prism.r2,workingDigits)/g.scale;
degree = vpa(pi,workingDigits)/180;
g.lat = vpa(point(1),workingDigits)*degree;
g.lon = vpa(point(2),workingDigits)*degree;
g.lat1 = vpa(prism.lat1,workingDigits)*degree;
g.lat2 = vpa(prism.lat2,workingDigits)*degree;
g.lon1 = vpa(prism.lon1,workingDigits)*degree;
g.lon2 = vpa(prism.lon2,workingDigits)*degree;
end

