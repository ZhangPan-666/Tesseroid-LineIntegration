function m = t1s_gauge_mean(a,r,branch)
switch branch
    case "source_outside_observer"
        m = 1-r^2/(3*a^2);
    case "source_inside_observer"
        m = 2*a/(3*r);
    otherwise
        error("Unsupported radial branch.");
end
end

