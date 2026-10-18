%% Statistical audit of the single-tesseroid runtime benchmark (Test No. 12)
% Reads the archived 160,801-point timing CSV, prints every statistic used
% in manuscript Section 4.1, and writes machine-readable summaries plus
% LaTeX tables for the Supplementary Information.

clc; clear;

here = fileparts(mfilename('fullpath'));
resultsDir = fullfile(here, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

inputFile = fullfile(here, 'RuntimeTest_LaplacePoisson_HeightPlus1km.csv');
logFile = fullfile(resultsDir, 'runtime_audit.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);

componentNames = ["V","Vx","Vy","Vz", ...
    "Vxx","Vxy","Vyy","Vzx","Vzy","Vzz", ...
    "Vxxx","Vxxy","Vxxz","Vxyz","Vyyx", ...
    "Vyyy","Vyyz","Vzzx","Vzzy","Vzzz"];
runtimeNames = componentNames + "_mean_s";

fprintf('============================================================\n');
fprintf('Test No. 12: single-tesseroid runtime statistical audit\n');
fprintf('============================================================\n');
fprintf('Input: %s\n', inputFile);

data = readtable(inputFile, 'VariableNamingRule', 'preserve');
required = ["longitude_deg","latitude_deg","inside_angular_bounds", ...
    runtimeNames,"all_components_mean_s"];
missing = setdiff(required, string(data.Properties.VariableNames));
assert(isempty(missing), 'Missing CSV variables: %s', strjoin(missing, ', '));

runtime = data{:, cellstr(runtimeNames)};
storedTotal = data.all_components_mean_s;
computedTotal = sum(runtime, 2);
inside = logical(data.inside_angular_bounds);
outside = ~inside;

assert(all(isfinite(runtime), 'all') && all(runtime > 0, 'all'), ...
    'All component times must be positive and finite.');
assert(all(isfinite(storedTotal)) && all(storedTotal > 0), ...
    'All stored cumulative times must be positive and finite.');
assert(any(inside) && any(outside), ...
    'Both inside- and outside-projection samples are required.');

lon = unique(data.longitude_deg);
lat = unique(data.latitude_deg);
assert(numel(lon)*numel(lat) == height(data), ...
    'The CSV does not contain a complete rectangular grid.');

fprintf('Grid: %d longitude x %d latitude = %d points\n', ...
    numel(lon), numel(lat), height(data));
fprintf('Inside projection: %d points; outside projection: %d points\n', ...
    nnz(inside), nnz(outside));
fprintf('Warm-up/timed repetitions: 2/10 (from benchmark source)\n');
fprintf('Maximum |stored total - component sum|: %.15e s\n\n', ...
    max(abs(storedTotal-computedTotal)));

statNames = {'Count','Minimum','Maximum','Mean','P50','P90','P99'};
componentTable = table('Size',[numel(componentNames), 8], ...
    'VariableTypes', [{'string'}, repmat({'double'},1,7)], ...
    'VariableNames', [{'Component'}, statNames]);

fprintf('Component statistics over all 160,801 points\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-6s %8s %12s %12s %12s %12s %12s %12s\n', ...
    'Name','N','Minimum','Maximum','Mean','P50','P90','P99');
for k = 1:numel(componentNames)
    s = stats(runtime(:,k));
    componentTable(k,:) = [{componentNames(k)}, num2cell(s)];
    fprintf('%-6s %8d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n', ...
        componentNames(k), s);
end

rankNames = ["potential","first-order sum","second-order sum", ...
    "third-order sum","all 20 components"];
rankValues = {runtime(:,1), sum(runtime(:,2:4),2), ...
    sum(runtime(:,5:10),2), sum(runtime(:,11:20),2), computedTotal};
regionNames = ["all points","inside projection","outside projection"];
regionMasks = {true(height(data),1), inside, outside};

rankTable = table('Size',[numel(rankNames)*numel(regionNames), 10], ...
    'VariableTypes', [{'string','string'}, repmat({'double'},1,8)], ...
    'VariableNames', {'Quantity','Region','Count','Minimum','Maximum','Mean', ...
    'P50','P90','P99','InsideOutsideMeanRatio'});

fprintf('\nDerivative-rank and cumulative statistics\n');
fprintf('------------------------------------------------------------\n');
row = 0;
for k = 1:numel(rankNames)
    insideMean = mean(rankValues{k}(inside));
    outsideMean = mean(rankValues{k}(outside));
    ratio = insideMean/outsideMean;
    for r = 1:numel(regionNames)
        row = row + 1;
        s = stats(rankValues{k}(regionMasks{r}));
        rankTable(row,:) = [{rankNames(k), regionNames(r)}, num2cell([s ratio])];
        fprintf('Test No. 12, %s, %s: N=%d, min=%.15e s, max=%.15e s, ', ...
            rankNames(k), regionNames(r), s(1), s(2), s(3));
        fprintf('mean=%.15e s, P50=%.15e s, P90=%.15e s, P99=%.15e s.\n', ...
            s(4), s(5), s(6), s(7));
    end
    fprintf('Test No. 12, %s: inside/outside mean ratio = %.8f.\n', ...
        rankNames(k), ratio);
end

storedStats = stats(storedTotal);
computedStats = stats(computedTotal);
[storedMax, storedMaxIndex] = max(storedTotal);
[componentMax, linearIndex] = max(runtime, [], 'all', 'linear');
[pointIndex, componentIndex] = ind2sub(size(runtime), linearIndex);

fprintf('\nAudit findings used by Section 4.1\n');
fprintf('------------------------------------------------------------\n');
fprintf('Stored cumulative time: min %.15e, max %.15e, mean %.15e, ', ...
    storedStats(2), storedStats(3), storedStats(4));
fprintf('P50 %.15e, P90 %.15e, P99 %.15e s.\n', ...
    storedStats(5), storedStats(6), storedStats(7));
fprintf('Recomputed cumulative time: min %.15e, max %.15e, mean %.15e, ', ...
    computedStats(2), computedStats(3), computedStats(4));
fprintf('P50 %.15e, P90 %.15e, P99 %.15e s.\n', ...
    computedStats(5), computedStats(6), computedStats(7));
fprintf('Largest stored cumulative time occurs at row %d: %.15e s.\n', ...
    storedMaxIndex, storedMax);
fprintf(['Maximum row geometry: longitude %.2f deg, latitude %.2f deg, ', ...
    'inside=%d, boundary distance %.3f deg.\n'], ...
    data.longitude_deg(storedMaxIndex), data.latitude_deg(storedMaxIndex), ...
    data.inside_angular_bounds(storedMaxIndex), ...
    data.distance_to_nearest_angular_boundary_deg(storedMaxIndex));
fprintf('Largest individual component time: %s at row %d = %.15e s.\n', ...
    componentNames(componentIndex), pointIndex, componentMax);
rankMeans = [mean(rankValues{1}), mean(rankValues{2}), ...
    mean(rankValues{3}), mean(rankValues{4})];
rankFractions = 100*rankMeans/mean(computedTotal);
fprintf(['Mean cumulative-time shares: potential %.3f%%, first order %.3f%%, ', ...
    'second order %.3f%%, third order %.3f%%.\n'], rankFractions);
fprintf('All component samples finite: yes.\n');

writetable(componentTable, fullfile(resultsDir, ...
    'runtime_component_statistics.csv'));
writetable(rankTable, fullfile(resultsDir, ...
    'runtime_rank_region_statistics.csv'));
writeLatex(fullfile(resultsDir, 'runtime_statistics_tables.tex'), ...
    componentTable, rankTable);

fprintf('\nGenerated outputs:\n');
fprintf('  %s\n', fullfile(resultsDir,'runtime_component_statistics.csv'));
fprintf('  %s\n', fullfile(resultsDir,'runtime_rank_region_statistics.csv'));
fprintf('  %s\n', fullfile(resultsDir,'runtime_statistics_tables.tex'));
fprintf('============================================================\n');
fprintf('Runtime audit completed successfully.\n');
fprintf('============================================================\n');
diary off;

function s = stats(x)
    x = x(:);
    s = [numel(x), min(x), max(x), mean(x), ...
        prctile(x,50), prctile(x,90), prctile(x,99)];
end

function writeLatex(fileName, componentTable, rankTable)
    fid = fopen(fileName, 'w');
    assert(fid ~= -1, 'Cannot write %s', fileName);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'%% Auto-generated by RunTimeTest1/Analyze_RuntimeTest1.m.\n');
    fprintf(fid,'\\subsection{Statistical audit of Test No.~12}\n');
    fprintf(fid,['The tables are regenerated from all 160,801 archived timing records. ', ...
        'P50, P90, and P99 are empirical percentiles; no timing sample is removed.\\par\n']);

    fprintf(fid,'\\begin{landscape}\n\\begin{table}[p]\\centering\\scriptsize\n');
    fprintf(fid,'\\setlength{\\tabcolsep}{4pt}\n');
    fprintf(fid,['\\caption{Per-component single-core timing statistics for Test No.~12. ', ...
        'Each point is the mean of ten timed evaluations after two warm-up evaluations.}', ...
        '\\label{tab:s-runtime-components}\n']);
    fprintf(fid,'\\resizebox{\\linewidth}{!}{%%\n');
    fprintf(fid,'\\begin{tabular}{lrrrrrrr}\\toprule\n');
    fprintf(fid,'Component & $N$ & Minimum & Maximum & Mean & P50 & P90 & P99 \\\\ \\midrule\n');
    for k=1:height(componentTable)
        fprintf(fid,'%s & %d & %s & %s & %s & %s & %s & %s \\\\\n', ...
            texComponent(componentTable.Component(k)), componentTable.Count(k), ...
            texSci(componentTable.Minimum(k)), texSci(componentTable.Maximum(k)), ...
            texSci(componentTable.Mean(k)), texSci(componentTable.P50(k)), ...
            texSci(componentTable.P90(k)), texSci(componentTable.P99(k)));
    end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n\\end{table}\n\\end{landscape}\n');

    allRows = rankTable(rankTable.Region=="all points",:);
    fprintf(fid,'\\begin{landscape}\n\\begin{table}[p]\\centering\\small\n');
    fprintf(fid,'\\setlength{\\tabcolsep}{4pt}\n');
    fprintf(fid,['\\caption{Derivative-rank and cumulative timing statistics for Test No.~12. ', ...
        'The final column compares mean costs inside and outside the tesseroid projection.}', ...
        '\\label{tab:s-runtime-ranks}\n']);
    fprintf(fid,'\\resizebox{\\linewidth}{!}{%%\n');
    fprintf(fid,'\\begin{tabular}{lrrrrrrrr}\\toprule\n');
    fprintf(fid,['Quantity & $N$ & Minimum & Maximum & Mean & P50 & P90 & P99 ', ...
        '& Inside/outside mean \\\\ \\midrule\n']);
    for k=1:height(allRows)
        fprintf(fid,'%s & %d & %s & %s & %s & %s & %s & %s & %.2f \\\\\n', ...
            allRows.Quantity(k), allRows.Count(k), texSci(allRows.Minimum(k)), ...
            texSci(allRows.Maximum(k)), texSci(allRows.Mean(k)), ...
            texSci(allRows.P50(k)), texSci(allRows.P90(k)), ...
            texSci(allRows.P99(k)), allRows.InsideOutsideMeanRatio(k));
    end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n\\end{table}\n\\end{landscape}\n');
    fprintf(fid,'\\FloatBarrier\n');
end

function out = texSci(x)
    if x==0, out='$0$'; return; end
    exponent = floor(log10(abs(x)));
    mantissa = x/10^exponent;
    out = sprintf('$%.2f\\times10^{%d}$',mantissa,exponent);
end

function out = texComponent(name)
    name = char(name);
    if strcmp(name,'V'), out='$V$'; else, out=sprintf('$V_{%s}$',name(2:end)); end
end
