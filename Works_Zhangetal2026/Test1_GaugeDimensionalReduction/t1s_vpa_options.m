function opts = t1s_vpa_options(workingDigits,quadratureDigits)
guarded = max(6,min(quadratureDigits,workingDigits-8));
opts.relTol = vpa("1e-"+string(guarded));
opts.absTol = opts.relTol;
end
