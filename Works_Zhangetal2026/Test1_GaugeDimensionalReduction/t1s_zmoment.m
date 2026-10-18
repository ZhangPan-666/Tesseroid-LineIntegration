function value = t1s_zmoment(r1,r2,s,c,n,p,q)
%T1S_ZMOMENT Radial-displacement-weighted normalized moment.

u1 = r1-c;
u2 = r2-c;
value = sym(0);
for k = 0:p
    coefR = nchoosek(p,k)*c^(p-k);
    for m = 0:q
        coefZ = nchoosek(q,m)*c^m*(-s^2)^(q-m);
        value = value+coefR*coefZ* ...
            t1s_base_moment(u1,u2,s,n,k+m);
    end
end
end

