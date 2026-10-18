function test_matlab(referenceCsv)
%TEST_MATLAB Check units/order against the separately built C++ SI example.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'Functions'));
model = [34.8 35.2 119.7 120.3 6341.2 6371.2 2670];
point = [35.11 120.08 6771.2;35.17 120.21 6471.2];
g = cell(1,20); gp=g;
[g{:}] = Tesseroid_GravityEstimate(model,point,1e-12,1e-9);
[gp{:}] = Tesseroid_GravityEstimateP(model,point,1e-12,1e-9);
assert(isequaln([g{:}],[gp{:}]),'Serial/OpenMP gravity mismatch');
twice=g;
[twice{:}] = Tesseroid_GravityEstimate([model;model],point,1e-12,1e-9);
assert(max(abs([twice{:}]-2*[g{:}]),[],'all') < 1e-10,'Superposition failed');
magModel = [model(1:6) 0.5 0.2 -0.1];
m=cell(1,10);mp=m;
[m{:}] = Tesseroid_MagneticEstimate(magModel,point,1e-12,1e-9);
[mp{:}] = Tesseroid_MagneticEstimateP(magModel,point,1e-12,1e-9);
assert(isequaln([m{:}],[mp{:}]),'Serial/OpenMP magnetic mismatch');
assert(all(isfinite([m{:}]),'all'),'Nonfinite magnetic output');
% One-based public-to-canonical permutation.
publicToCanonical = [1 2 3 4 5 6 8 7 9 10 11 12 13 15 14 17 18 16 19 20];
orders = [0 ones(1,3) 2*ones(1,6) 3*ones(1,10)];
if nargin > 0
    reference = readtable(referenceCsv,'VariableNamingRule','preserve');
    for j=1:20
        k=publicToCanonical(j);
        si=g{j}(1)*6.67430e-11*1000^(2-orders(k));
        assert(abs(si-reference.value(k)) < 2e-7*max(abs(reference.value(k)),1e-25), ...
            'C++/MEX unit or component mismatch at output %d.',j);
    end
end
rejected = false;
try
    [g{:}] = Tesseroid_GravityEstimate(model,point,-1,1e-9);
catch errorInfo
    rejected=strcmp(errorInfo.identifier,'TesseroidLineIntegration:InvalidTolerance');
end
assert(rejected,'Negative tolerance was not rejected');
fprintf('MATLAB gravity/magnetic serial/OpenMP, superposition, order/units and input checks passed.\n');
end
