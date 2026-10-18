%% Test 3: polar analytical validation for four tesseroids
% The reference is Tesseroid_AnalyticalPolar_vpa.m.  The numerical result
% is the ordinary Tesseroid gravity estimator evaluated exactly at the
% North Pole, where the polar branch is exercised.

clc; clear; close all;

here = fileparts(mfilename('fullpath'));
funcDir = getenv('TLI_POLAR_REFERENCE_DIR');
assert(isfolder(funcDir), 'Set TLI_POLAR_REFERENCE_DIR to the external analytical-reference Functions directory.');
addpath(funcDir);
assert(exist('Tesseroid_AnalyticalPolar_vpa','file') == 2, 'Missing external Tesseroid_AnalyticalPolar_vpa.m.');

RE = 6371.2;                         % km
height = (-1000:1:1000).';           % km relative to RE
radius = RE + height;
lambdaObs = zeros(size(radius));     % longitude is immaterial at the pole
latitudeObs = 90 * ones(size(radius));
precisionDigits = 128;
absoluteTolerance = 1e-12;
relativeTolerance = 1e-12;

% [Fai1 Fai2 Lamda1 Lamda2 R1 R2] in degrees and km.
cases = struct(...
    'name', {'Normal','Small','Thin','ZonalBand'}, ...
    'label', {'(a)','(b)','(c)','(d)'}, ...
    'parameters', {...
        [69.99, 89.99, -10, 10, 6341.2, 6371.2], ...
        [89.89, 89.99, -0.05, 0.05, 6370.2, 6371.2], ...
        [69.99, 89.99, -10, 10, 6371.19, 6371.2], ...
        [69.99, 89.99, -180, 180, 6341.2, 6371.2]});

componentNames = {'V','Vz','Vxx','Vyy','Vzz','Vxxz','Vyyz','Vzzz'};
numCase = numel(cases);
numPoint = numel(height);
numComponent = numel(componentNames);
analytical = zeros(numPoint, numComponent, numCase);
cppTable = readtable(fullfile(here,'results','polar_cpp_results.csv'), ...
    'VariableNamingRule','preserve');
numerical = zeros(numPoint, numComponent, numCase);
absoluteError = zeros(numPoint, numComponent, numCase);
relativeError = zeros(numPoint, numComponent, numCase);

for k = 1:numCase
    p = cases(k).parameters;
    fprintf('Case %s: %s\n', cases(k).label, cases(k).name);

    [Va,Vza,Vxxa,Vyya,Vzza,Vxxza,Vyyza,Vzzza] = ...
        Tesseroid_AnalyticalPolar_vpa(p(1),p(2),p(3),p(4), ...
        lambdaObs,p(5),p(6),radius,precisionDigits);
    reference = double([Va,Vza,Vxxa,Vyya,Vzza,Vxxza,Vyyza,Vzzza]);

    rows = find(strcmp(string(cppTable.('case')), cases(k).name));
    if numel(rows) ~= numPoint || any(cppTable.height(rows) ~= height)
        error('C++ result does not contain the expected 1-km profile for %s.', cases(k).name);
    end
    estimate = [cppTable.V(rows),cppTable.Vz(rows),cppTable.Vxx(rows),cppTable.Vyy(rows), ...
        cppTable.Vzz(rows),cppTable.Vxxz(rows),cppTable.Vyyz(rows),cppTable.Vzzz(rows)];

    analytical(:,:,k) = reference;
    numerical(:,:,k) = estimate;
    absoluteError(:,:,k) = abs(estimate-reference);
    scale = max(abs(reference), realmin('double'));
    relativeError(:,:,k) = absoluteError(:,:,k)./scale;

    cases(k).maxAbsoluteError = max(absoluteError(:,:,k),[],1);
    cases(k).rmseAbsoluteError = sqrt(mean(absoluteError(:,:,k).^2,1));
    cases(k).maxRelativeError = max(relativeError(:,:,k),[],1);
    cases(k).rmseRelativeError = sqrt(mean(relativeError(:,:,k).^2,1));
end

% Write a compact, machine-readable summary.
summaryCase = strings(numCase*numComponent,1);
summaryComponent = strings(numCase*numComponent,1);
maxAbs = zeros(size(summaryCase)); rmseAbs = zeros(size(summaryCase));
maxRel = zeros(size(summaryCase)); rmseRel = zeros(size(summaryCase));
q = 0;
for k = 1:numCase
    for j = 1:numComponent
        q = q + 1;
        summaryCase(q) = cases(k).name;
        summaryComponent(q) = componentNames{j};
        maxAbs(q) = cases(k).maxAbsoluteError(j);
        rmseAbs(q) = cases(k).rmseAbsoluteError(j);
        maxRel(q) = cases(k).maxRelativeError(j);
        rmseRel(q) = cases(k).rmseRelativeError(j);
    end
end
summary = table(summaryCase,summaryComponent,maxAbs,rmseAbs,maxRel,rmseRel, ...
    'VariableNames',{'case','component','max_abs','rmse_abs','max_rel','rmse_rel'});
writetable(summary,fullfile(here,'results','polar_summary.csv'));
save(fullfile(here,'results','polar_validation.mat'), ...
    'cases','height','radius','componentNames','analytical','numerical', ...
    'absoluteError','relativeError','summary','-v7.3');

fprintf('\nWrote results/polar_summary.csv and results/polar_validation.mat\n');
disp(summary);
