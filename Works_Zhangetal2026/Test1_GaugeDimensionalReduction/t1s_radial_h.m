function h = t1s_radial_h(r1,r2,r,c)
s = r*sqrt(1-c^2);
u1 = r1-r*c;
u2 = r2-r*c;
l1 = sqrt(u1^2+s^2);
l2 = sqrt(u2^2+s^2);
h = r*((l2-l1)+r*c*(asinh(u2/s)-asinh(u1/s)));
end

