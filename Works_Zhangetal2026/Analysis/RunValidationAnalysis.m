%% Unified statistical audit for manuscript Validation Tests No. 1--11
% This script reads the archived numerical outputs, recomputes every
% statistic quoted in the Validation section, prints an auditable report,
% and writes machine-readable CSV files plus LaTeX tables for the SI.

clc; clear;

here = fileparts(mfilename('fullpath'));
worksDir = fileparts(here);
resultsDir = fullfile(here, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

logFile = fullfile(resultsDir, 'validation_audit.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);
fprintf('Unified validation audit for Tests No. 1--11\n');
fprintf('Generated: %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('Statistics: minimum, maximum, arithmetic mean, RMSE, P50, P90, P99.\n\n');

names = {'TestNo','Family','Configuration','Component','ErrorType', ...
    'Count','NaNCount','Minimum','Maximum','Mean','RMSE','P50','P90','P99'};
types = {'double','string','string','string','string', ...
    'double','double','double','double','double','double','double','double','double'};
detail = table('Size',[0 numel(names)],'VariableTypes',types,'VariableNames',names);
summary = detail;
claims = table('Size',[0 3], ...
    'VariableTypes',{'string','double','string'}, ...
    'VariableNames',{'Key','Value','Description'});

%% Test No. 1: dimensional reduction
fprintf('============================================================\n');
fprintf('Test No. 1: direct dimensional-reduction verification\n');
fprintf('============================================================\n');
t1Root = fullfile(worksDir, 'Test1_GaugeDimensionalReduction', 'results');
heights = [100 400];
comparisons = ["2-D vs 3-D", "1-D vs 3-D", "double-precision TLI vs 3-D"];
variables = ["rel3d2d", "rel3d1d", "rel3dCpp"];
allT1 = cell(numel(heights),1);
for ih = 1:numel(heights)
    T = readtable(fullfile(t1Root, sprintf('plus%dkm',heights(ih)), ...
        'summary_current.csv'), 'VariableNamingRule','preserve');
    assert(height(T)==20, 'Test No. 1 must contain 20 components.');
    assert(all(T.done1d==1 & T.done2d==1 & T.done3d==1), ...
        'A variable-precision path is incomplete.');
    assert(all(T.conv1d==1 & T.conv2d==1 & T.conv3d==1), ...
        'A variable-precision path did not converge.');
    allT1{ih} = T;
    for ip = 1:numel(comparisons)
        x = T.(char(variables(ip)));
        detail = addRow(detail,1,"Dimensional reduction", ...
            sprintf('%d km',heights(ih)),comparisons(ip),"relative difference",x);
        printStats(sprintf('Test No. 1, %d km, %s',heights(ih),comparisons(ip)),x);
    end
    t1d = sum(T.seconds1d); t2d = sum(T.seconds2d); t3d = sum(T.seconds3d);
    speed = t3d/t1d;
    fprintf('Test No. 1, %d km: summed time 1-D=%.1f s, 2-D=%.1f s, 3-D=%.1f s; 3-D/1-D=%.2f.\n', ...
        heights(ih),t1d,t2d,t3d,speed);
    claims = addClaim(claims,sprintf('test1_%dkm_time_1d',heights(ih)),t1d,'Summed 1-D variable-precision runtime (s)');
    claims = addClaim(claims,sprintf('test1_%dkm_time_2d',heights(ih)),t2d,'Summed 2-D variable-precision runtime (s)');
    claims = addClaim(claims,sprintf('test1_%dkm_time_3d',heights(ih)),t3d,'Summed 3-D variable-precision runtime (s)');
    claims = addClaim(claims,sprintf('test1_%dkm_speedup',heights(ih)),speed,'3-D runtime divided by 1-D runtime');
end
for ip = 1:numel(comparisons)
    x = [allT1{1}.(char(variables(ip))); allT1{2}.(char(variables(ip)))];
    s = getStats(x);
    claims = addClaim(claims,"test1_overall_max_"+variables(ip),s.Maximum, ...
        "Maximum componentwise relative difference across both heights: "+comparisons(ip));
end

%% Tests No. 2--5: polar-axis analytical validation
fprintf('\n============================================================\n');
fprintf('Tests No. 2--5: analytical polar-axis validation\n');
fprintf('============================================================\n');
polarFile = fullfile(worksDir,'Test3_PolarValidation','results','polar_validation.mat');
P = load(polarFile,'absoluteError','relativeError','componentNames','cases','height');
assert(isequal(size(P.absoluteError),[2001 8 4]), 'Unexpected polar absolute-error array size.');
assert(isequal(size(P.relativeError),[2001 8 4]), 'Unexpected polar relative-error array size.');
for k = 1:4
    testNo = k+1;
    caseName = string(P.cases(k).name);
    rmseRel = zeros(8,1);
    for j = 1:8
        comp = string(P.componentNames{j});
        xa = P.absoluteError(:,j,k);
        xr = P.relativeError(:,j,k);
        detail = addRow(detail,testNo,"Polar axis",caseName,comp,"absolute error",xa);
        detail = addRow(detail,testNo,"Polar axis",caseName,comp,"relative error",xr);
        printStats(sprintf('Test No. %d (%s), %s absolute error',testNo,caseName,comp),xa);
        printStats(sprintf('Test No. %d (%s), %s relative error',testNo,caseName,comp),xr);
        rmseRel(j) = getStats(xr).RMSE;
    end
    [worst,idx] = max(rmseRel);
    fprintf('Test No. %d (%s): largest componentwise relative-error RMSE = %.15e (%s).\n', ...
        testNo,caseName,worst,P.componentNames{idx});
    claims = addClaim(claims,sprintf('test%d_worst_component_rel_rmse',testNo),worst, ...
        "Largest componentwise relative-error RMSE; component="+string(P.componentNames{idx}));
    pooledAbs = reshape(P.absoluteError(:,:,k),[],1);
    pooledRel = reshape(P.relativeError(:,:,k),[],1);
    summary = addRow(summary,testNo,"Polar axis",caseName,"all eight pooled","absolute error",pooledAbs);
    summary = addRow(summary,testNo,"Polar axis",caseName,"all eight pooled","relative error",pooledRel);
end

%% Tests No. 6--8: complete spherical shell
fprintf('\n============================================================\n');
fprintf('Tests No. 6--8: complete-spherical-shell validation\n');
fprintf('============================================================\n');
shellRoot = fullfile(worksDir,'Test2_SphericalShellValidation','results');
shellFiles = ["ResultsPan_10x10Split_1km.dat", ...
              "ResultsPan_10x10Split_-31km.dat", ...
              "ResultsPan_10x10Split_-15km.dat"];
shellLabels = ["outside, R_E+1 km", "cavity, R_E-31 km", "material, R_E-15 km"];
componentNames = ["V","Vx","Vy","Vz","Vxx","Vxy","Vxz","Vyy","Vyz","Vzz", ...
    "Vxxx","Vxxy","Vxxz","Vxyy","Vxyz","Vxzz","Vyyy","Vyyz","Vyzz","Vzzz"];
canonicalFromWork = [1 2 3 4 5 6 8 7 9 10 11 12 13 15 14 18 16 17 19 20];
orders = [0 ones(1,3) 2*ones(1,6) 3*ones(1,10)];
for k = 1:3
    testNo = k+5;
    [lon,lat,r,workValues] = readWorkFile(fullfile(shellRoot,shellFiles(k)));
    numerical = workValues(:,canonicalFromWork);
    reference = shellReference(r(1));
    errors = abs(numerical-reference);
    assert(numel(lon)==64800 && numel(lat)==64800, 'Shell test must contain 64,800 points.');
    for j = 1:20
        x = errors(:,j);
        detail = addRow(detail,testNo,"Spherical shell",shellLabels(k), ...
            componentNames(j),"absolute error",x);
        printStats(sprintf('Test No. %d (%s), %s absolute error',testNo,shellLabels(k),componentNames(j)),x);
    end
    for ord = 0:3
        pooled = reshape(errors(:,orders==ord),[],1);
        summary = addRow(summary,testNo,"Spherical shell",shellLabels(k), ...
            sprintf('derivative order %d pooled',ord),"absolute error",pooled);
    end
    nanCount = sum(~isfinite(errors),'all');
    fprintf('Test No. %d: NaN/Inf error count = %d of %d comparisons.\n',testNo,nanCount,numel(errors));
    claims = addClaim(claims,sprintf('test%d_nonfinite_count',testNo),nanCount,'Nonfinite shell-error count');
    for ord = 0:3
        compMax = max(errors(:,orders==ord),[],1);
        [rankMax,ii] = max(compMax);
        rankNames = componentNames(orders==ord);
        claims = addClaim(claims,sprintf('test%d_order%d_max_abs',testNo,ord),rankMax, ...
            "Largest pointwise absolute error at derivative order "+ord+"; component="+rankNames(ii));
        fprintf('Test No. %d: maximum order-%d absolute error = %.15e (%s).\n', ...
            testNo,ord,rankMax,rankNames(ii));
    end
    if testNo==7
        derivativeExactZero = all(numerical(:,2:end)==0,'all');
        claims = addClaim(claims,'test7_all_derivatives_exact_zero',double(derivativeExactZero), ...
            'One if all 19 derivative outputs are exactly zero in double precision');
        claims = addClaim(claims,'test7_potential_abs_error',errors(1,1), ...
            'Absolute potential error in the shell cavity');
        fprintf('Test No. 7: all 19 derivative outputs exact zero = %d; potential error = %.15e.\n', ...
            derivativeExactZero,errors(1,1));
    end
end

%% Tests No. 9--11: Laplace--Poisson and differentiated traces
fprintf('\n============================================================\n');
fprintf('Tests No. 9--11: Laplace--Poisson identity residuals\n');
fprintf('============================================================\n');
traceRoot = fullfile(worksDir,'Test4_LaplacePoissonValidation');
traceFiles = ["ResultsPan_1km.dat","ResultsPan_-31km.dat","ResultsPan_-15km.dat"];
traceLabels = ["above, R_E+1 km","below, R_E-31 km","radially inside, R_E-15 km"];
identityNames = ["Vxx+Vyy+Vzz", "Vxxx+Vxyy+Vxzz", ...
                 "Vxxy+Vyyy+Vyzz", "Vxxz+Vyyz+Vzzz"];
for k = 1:3
    testNo = k+8;
    [lon,lat,~,F] = readWorkFile(fullfile(traceRoot,traceFiles(k)));
    inside = lon>110 & lon<130 & lat>30 & lat<50;
    theory2 = zeros(size(lon));
    if testNo==11, theory2(inside) = -4*pi; end
    residuals = zeros(numel(lon),4);
    residuals(:,1) = abs(F(:,5)+F(:,7)+F(:,10)-theory2);
    residuals(:,2) = abs(F(:,11)+F(:,15)+F(:,18));
    residuals(:,3) = abs(F(:,12)+F(:,16)+F(:,19));
    residuals(:,4) = abs(F(:,13)+F(:,17)+F(:,20));
    assert(numel(lon)==160801, 'Trace test must contain 160,801 points.');
    for j = 1:4
        x = residuals(:,j);
        detail = addRow(detail,testNo,"Trace identities",traceLabels(k), ...
            identityNames(j),"absolute residual",x);
        printStats(sprintf('Test No. %d (%s), %s residual',testNo,traceLabels(k),identityNames(j)),x);
    end
    summary = addRow(summary,testNo,"Trace identities",traceLabels(k), ...
        "all four pooled","absolute residual",residuals(:));
    secondMax = max(residuals(:,1));
    [thirdMax,linearIdx] = max(residuals(:,2:4),[],'all','linear');
    [~,thirdIdentity] = ind2sub(size(residuals(:,2:4)),linearIdx);
    countAbove = sum(residuals(:,1)>1e-8);
    claims = addClaim(claims,sprintf('test%d_second_trace_max',testNo),secondMax, ...
        'Maximum absolute Hessian-trace residual');
    claims = addClaim(claims,sprintf('test%d_largest_differentiated_trace_max',testNo),thirdMax, ...
        "Largest differentiated-trace residual; identity="+identityNames(thirdIdentity+1));
    claims = addClaim(claims,sprintf('test%d_second_trace_count_gt_1e8',testNo),countAbove, ...
        'Number of Hessian-trace residuals greater than 1e-8');
    fprintf('Test No. %d: max Hessian-trace residual=%.15e; count >1e-8=%d; largest differentiated-trace residual=%.15e (%s).\n', ...
        testNo,secondMax,countAbove,thirdMax,identityNames(thirdIdentity+1));
end

%% Write outputs
writetable(detail,fullfile(resultsDir,'validation_component_statistics.csv'));
writetable(summary,fullfile(resultsDir,'validation_experiment_statistics.csv'));
writetable(claims,fullfile(resultsDir,'manuscript_claims.csv'));
writeLatexTables(detail,fullfile(resultsDir,'validation_statistics_tables.tex'));

fprintf('\n============================================================\n');
fprintf('Audit complete.\n');
fprintf('Detailed statistics: %s\n',fullfile(resultsDir,'validation_component_statistics.csv'));
fprintf('Experiment summaries: %s\n',fullfile(resultsDir,'validation_experiment_statistics.csv'));
fprintf('Manuscript claims: %s\n',fullfile(resultsDir,'manuscript_claims.csv'));
fprintf('SI tables: %s\n',fullfile(resultsDir,'validation_statistics_tables.tex'));
fprintf('============================================================\n');
diary off;

%% Local functions
function T = addRow(T,testNo,family,configuration,component,errorType,x)
    s = getStats(x);
    T(end+1,:) = {testNo,string(family),string(configuration),string(component), ...
        string(errorType),s.Count,s.NaNCount,s.Minimum,s.Maximum,s.Mean,s.RMSE, ...
        s.P50,s.P90,s.P99};
end

function C = addClaim(C,key,value,description)
    C(end+1,:) = {string(key),double(value),string(description)};
end

function s = getStats(x)
    x = double(x(:));
    s.NaNCount = sum(~isfinite(x));
    x = x(isfinite(x));
    s.Count = numel(x);
    if isempty(x)
        [s.Minimum,s.Maximum,s.Mean,s.RMSE,s.P50,s.P90,s.P99] = deal(NaN);
        return;
    end
    s.Minimum = min(x);
    s.Maximum = max(x);
    s.Mean = mean(x);
    s.RMSE = norm(x)/sqrt(numel(x));
    q = prctile(x,[50 90 99]);
    s.P50 = q(1); s.P90 = q(2); s.P99 = q(3);
end

function printStats(label,x)
    s = getStats(x);
    fprintf('%s: N=%d, NaN/Inf=%d, min=%.15e, max=%.15e, mean=%.15e, RMSE=%.15e, P50=%.15e, P90=%.15e, P99=%.15e.\n', ...
        label,s.Count,s.NaNCount,s.Minimum,s.Maximum,s.Mean,s.RMSE,s.P50,s.P90,s.P99);
end

function [lon,lat,r,F] = readWorkFile(path)
    fid = fopen(path,'rb');
    assert(fid>=0,'Cannot open %s.',path);
    cleaner = onCleanup(@()fclose(fid));
    n = fread(fid,1,'int32=>double');
    assert(~isempty(n) && n>0,'Invalid Work file header: %s.',path);
    lon = fread(fid,n,'double=>double');
    lat = fread(fid,n,'double=>double');
    r = fread(fid,n,'double=>double');
    F = zeros(n,20);
    for j = 1:20, F(:,j) = fread(fid,n,'double=>double'); end
    assert(all(size(F)==[n 20]) && numel(lon)==n && numel(lat)==n && numel(r)==n, ...
        'Truncated Work file: %s.',path);
    extra = fread(fid,1,'uint8');
    assert(isempty(extra),'Unexpected trailing bytes in %s.',path);
end

function a = shellReference(r)
    r1 = 6341.2; r2 = 6371.2;
    d2 = (r2-r1)*(r2+r1);
    d3 = (r2-r1)*(r2^2+r2*r1+r1^2);
    a = zeros(1,20);
    if r >= r2
        V = 4*pi*d3/(3*r); Vz = -4*pi*d3/(3*r^2);
        Vxx = -4*pi*d3/(3*r^3); Vzz = 8*pi*d3/(3*r^3);
        Vxxz = 4*pi*d3/r^4; Vzzz = -8*pi*d3/r^4;
    elseif r <= r1
        V = 2*pi*d2; Vz = 0; Vxx = 0; Vzz = 0; Vxxz = 0; Vzzz = 0;
    else
        dr = (r-r1)*(r^2+r*r1+r1^2);
        V = 2*pi*(d2-(r-r1)^2*(r+2*r1)/(3*r));
        Vz = -4*pi*dr/(3*r^2); Vxx = -4*pi*dr/(3*r^3);
        Vzz = -4*pi/3*(1+2*r1^3/r^3);
        Vxxz = -4*pi*r1^3/r^4; Vzzz = 8*pi*r1^3/r^4;
    end
    a([1 4 5 8 10 13 18 20]) = [V Vz Vxx Vxx Vzz Vxxz Vxxz Vzzz];
end

function writeLatexTables(T,path)
    fid = fopen(path,'w');
    assert(fid>=0,'Cannot write %s.',path);
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid,'%% Auto-generated by Analysis/RunValidationAnalysis.m. Do not edit manually.\n');
    fprintf(fid,'\\subsection{Reproducible statistical audit of Tests No.~1--11}\n');
    fprintf(fid,['All entries below are recomputed from the archived pointwise outputs. ' ...
        'For a nonnegative error sample $e_i$, RMSE is $[N^{-1}\\sum_i e_i^2]^{1/2}$; ' ...
        'P50, P90, and P99 are empirical percentiles. No finite value is removed. ' ...
        'The accompanying MATLAB audit reports the count of nonfinite values separately.\\par\n']);

    rows = T(T.TestNo==1,:);
    tableBegin(fid,'Statistics of the componentwise relative differences in Test No.~1.', ...
        'tab:s-audit-test1','llrrrrrrrr');
    fprintf(fid,'Height & Comparison & $N$ & Minimum & Maximum & Mean & RMSE & P50 & P90 & P99 \\\\ \\midrule\n');
    for i=1:height(rows)
        fprintf(fid,'%s & %s & %d & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
            rows.Configuration(i),rows.Component(i),rows.Count(i), ...
            texSci(rows.Minimum(i)),texSci(rows.Maximum(i)),texSci(rows.Mean(i)), ...
            texSci(rows.RMSE(i)),texSci(rows.P50(i)),texSci(rows.P90(i)),texSci(rows.P99(i)));
    end
    tableEnd(fid);

    for testNo=2:5
        rows = T(T.TestNo==testNo,:);
        caption = sprintf('Componentwise absolute- and relative-error statistics for Test No.~%d over 2,001 polar-axis points.',testNo);
        label = sprintf('tab:s-audit-polar%d',testNo);
        tableBegin(fid,caption,label,'llrrrrrrrr');
        fprintf(fid,'Component & Error & $N$ & Minimum & Maximum & Mean & RMSE & P50 & P90 & P99 \\\\ \\midrule\n');
        for i=1:height(rows)
            fprintf(fid,'%s & %s & %d & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
                texComponent(rows.Component(i)),shortError(rows.ErrorType(i)),rows.Count(i), ...
                texSci(rows.Minimum(i)),texSci(rows.Maximum(i)),texSci(rows.Mean(i)), ...
                texSci(rows.RMSE(i)),texSci(rows.P50(i)),texSci(rows.P90(i)),texSci(rows.P99(i)));
        end
        tableEnd(fid);
    end

    for testNo=6:8
        rows = T(T.TestNo==testNo,:);
        caption = sprintf('Componentwise absolute-error statistics for the complete-spherical-shell Test No.~%d over 64,800 points.',testNo);
        label = sprintf('tab:s-audit-shell%d',testNo);
        tableBegin(fid,caption,label,'lrrrrrrrr');
        fprintf(fid,'Component & $N$ & Minimum & Maximum & Mean & RMSE & P50 & P90 & P99 \\\\ \\midrule\n');
        for i=1:height(rows)
            fprintf(fid,'%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
                texComponent(rows.Component(i)),rows.Count(i),texSci(rows.Minimum(i)), ...
                texSci(rows.Maximum(i)),texSci(rows.Mean(i)),texSci(rows.RMSE(i)), ...
                texSci(rows.P50(i)),texSci(rows.P90(i)),texSci(rows.P99(i)));
        end
        tableEnd(fid);
    end

    for testNo=9:11
        rows = T(T.TestNo==testNo,:);
        caption = sprintf('Absolute-residual statistics for the Laplace--Poisson identities in Test No.~%d over 160,801 points.',testNo);
        label = sprintf('tab:s-audit-trace%d',testNo);
        tableBegin(fid,caption,label,'lrrrrrrrr');
        fprintf(fid,'Identity & $N$ & Minimum & Maximum & Mean & RMSE & P50 & P90 & P99 \\\\ \\midrule\n');
        for i=1:height(rows)
            fprintf(fid,'%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
                texIdentity(rows.Component(i)),rows.Count(i),texSci(rows.Minimum(i)), ...
                texSci(rows.Maximum(i)),texSci(rows.Mean(i)),texSci(rows.RMSE(i)), ...
                texSci(rows.P50(i)),texSci(rows.P90(i)),texSci(rows.P99(i)));
        end
        tableEnd(fid);
    end
    fprintf(fid,'\\FloatBarrier\n');
end

function tableBegin(fid,caption,label,columns)
    fprintf(fid,'\\begin{landscape}\n\\begin{table}[p]\\centering\\scriptsize\n');
    fprintf(fid,'\\setlength{\\tabcolsep}{3.2pt}\n');
    fprintf(fid,'\\caption{%s}\\label{%s}\n',caption,label);
    fprintf(fid,'\\resizebox{\\linewidth}{!}{%%\n\\begin{tabular}{%s}\n\\toprule\n',columns);
end

function tableEnd(fid)
    fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table}\n\\end{landscape}\n');
end

function out = texSci(x)
    if isnan(x), out='--'; return; end
    if x==0, out='$0$'; return; end
    exponent = floor(log10(abs(x)));
    mantissa = x/10^exponent;
    out = sprintf('$%.2f\\times10^{%d}$',mantissa,exponent);
end

function out = texComponent(name)
    name = char(name);
    if strcmp(name,'V'), out='$V$'; else, out=sprintf('$V_{%s}$',name(2:end)); end
end

function out = texIdentity(name)
    parts = split(string(name),'+');
    terms = strings(size(parts));
    for i=1:numel(parts), terms(i)=string(texComponent(parts(i))); end
    out = char(join(terms,'+'));
    out = strrep(out,'$+$','+');
end

function out = shortError(name)
    if startsWith(name,'absolute'), out='absolute'; else, out='relative'; end
end
