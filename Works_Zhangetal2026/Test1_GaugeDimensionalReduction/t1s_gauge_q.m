function q = t1s_gauge_q(a,r,c,branch)
m = t1s_gauge_mean(a,r,branch);
d = sqrt((a-r)^2+2*a*r*(1-c));
p = -((a^2-r^2)*d+d^3/3)/(2*a^2*r);
dm = a+r;
pm = -((a^2-r^2)*dm+dm^3/3)/(2*a^2*r);
q = a^2/(1-c^2)*(p-pm-m*(c+1));
end

