function value = t1s_moment(r1,r2,s,c,n,p)
%T1S_MOMENT Ordinary normalized radial moment.

u1 = r1-c;
u2 = r2-c;
value = sym(0);
for k = 0:p
    value = value+nchoosek(p,k)*c^(p-k)* ...
        t1s_base_moment(u1,u2,s,n,k);
end
end

