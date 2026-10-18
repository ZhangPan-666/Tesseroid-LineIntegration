function [tx,ty,cosine] = t1s_geometry_terms(sourceLat,sourceLon,g)
%T1S_GEOMETRY_TERMS Source direction in the observation local frame.

beta = sourceLon-g.lon;
tx = cos(g.lat)*sin(sourceLat)- ...
    sin(g.lat)*cos(sourceLat)*cos(beta);
ty = cos(sourceLat)*sin(beta);
cosine = sin(g.lat)*sin(sourceLat)+ ...
    cos(g.lat)*cos(sourceLat)*cos(beta);
end

