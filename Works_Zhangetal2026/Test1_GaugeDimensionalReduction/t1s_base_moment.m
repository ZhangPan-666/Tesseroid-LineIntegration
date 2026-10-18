function b = t1s_base_moment(u1,u2,s,n,k)
%T1S_BASE_MOMENT Exact endpoint recurrence for integral u^k/(u^2+s^2)^(n/2).

b = local_primitive(u2,s,n,k)-local_primitive(u1,s,n,k);
end

function f = local_primitive(u,s,n,k)
if k >= 2
    f = local_primitive(u,s,n-2,k-2)-s^2*local_primitive(u,s,n,k-2);
    return
end

l = sqrt(u^2+s^2);
if k == 1
    if n == 1
        f = l;
    else
        f = l^(2-n)/(2-n);
    end
    return
end

if n == 1
    f = asinh(u/s);
elseif n < 1
    % Here -n is a positive odd power of l=sqrt(u^2+s^2).
    power = -n;
    f = u*l^power/(power+1)+power*s^2/(power+1)* ...
        local_primitive(u,s,n+2,0);
else
    f = u/((n-2)*s^2*l^(n-2))+ ...
        (n-3)/((n-2)*s^2)*local_primitive(u,s,n-2,0);
end
end
