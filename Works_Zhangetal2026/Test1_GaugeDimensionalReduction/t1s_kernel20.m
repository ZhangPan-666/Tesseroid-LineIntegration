function kernels = t1s_kernel20(dx,dy,dz,distance,jacobian)
%T1S_KERNEL20 Cartesian observation derivatives of the Newton kernel.

d = [dx,dy,dz];
kernels = sym(zeros(20,1));
kernels(1) = jacobian/distance;
for i = 1:3
    kernels(1+i) = jacobian*d(i)/distance^3;
end

secondPairs = [1 1;1 2;1 3;2 2;2 3;3 3];
for p = 1:6
    i = secondPairs(p,1);
    j = secondPairs(p,2);
    delta = double(i==j);
    kernels(4+p) = jacobian*( ...
        3*d(i)*d(j)/distance^5-delta/distance^3);
end

thirdTriples = [1 1 1;1 1 2;1 1 3;1 2 2;1 2 3;1 3 3; ...
                2 2 2;2 2 3;2 3 3;3 3 3];
for p = 1:10
    i = thirdTriples(p,1);
    j = thirdTriples(p,2);
    k = thirdTriples(p,3);
    deltaTerm = double(i==j)*d(k)+double(i==k)*d(j)+ ...
        double(j==k)*d(i);
    kernels(10+p) = jacobian*( ...
        15*d(i)*d(j)*d(k)/distance^7-3*deltaTerm/distance^5);
end
end

