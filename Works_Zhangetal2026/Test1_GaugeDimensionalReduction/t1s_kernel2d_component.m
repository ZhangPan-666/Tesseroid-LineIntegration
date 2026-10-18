function kernel = t1s_kernel2d_component(g,phi,lambda,index)
%T1S_KERNEL2D_COMPONENT Lazily construct one radial-moment kernel.

[tx,ty,c] = t1s_geometry_terms(phi,lambda,g);
s = sqrt(tx^2+ty^2);
r1 = g.r1;
r2 = g.r2;
M = @(n,p)t1s_moment(r1,r2,s,c,n,p);
Z = @(n,p,q)t1s_zmoment(r1,r2,s,c,n,p,q);

switch index
    case 1,  q=M(1,2);
    case 2,  q=tx*M(3,3);
    case 3,  q=ty*M(3,3);
    case 4,  q=Z(3,2,1);
    case 5,  q=3*tx^2*M(5,4)-M(3,2);
    case 6,  q=3*tx*ty*M(5,4);
    case 7,  q=3*tx*Z(5,3,1);
    case 8,  q=3*ty^2*M(5,4)-M(3,2);
    case 9,  q=3*ty*Z(5,3,1);
    case 10, q=3*Z(5,2,2)-M(3,2);
    case 11, q=15*tx^3*M(7,5)-9*tx*M(5,3);
    case 12, q=15*tx^2*ty*M(7,5)-3*ty*M(5,3);
    case 13, q=15*tx^2*Z(7,4,1)-3*Z(5,2,1);
    case 14, q=15*ty^2*tx*M(7,5)-3*tx*M(5,3);
    case 15, q=15*tx*ty*Z(7,4,1);
    case 16, q=15*tx*Z(7,3,2)-3*tx*M(5,3);
    case 17, q=15*ty^3*M(7,5)-9*ty*M(5,3);
    case 18, q=15*ty^2*Z(7,4,1)-3*Z(5,2,1);
    case 19, q=15*ty*Z(7,3,2)-3*ty*M(5,3);
    case 20, q=15*Z(7,2,3)-9*Z(5,2,1);
    otherwise, error("Component index must be from 1 to 20.");
end
kernel = cos(phi)*q;
end

