% Self-contained scaled-unit example. Run from any directory after building MEX.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'Functions'));
model = [34.8 35.2 119.7 120.3 6341.2 6371.2 2670];
points = [35.11 120.08 6771.2; 35.17 120.21 6471.2];
gravity = cell(1,20);
[gravity{:}] = Tesseroid_GravityEstimate(model,points,1e-12,1e-9);
disp(table(points(:,1),points(:,2),gravity{1},gravity{4}, ...
    'VariableNames',{'latitude','longitude','V_scaled','Vz_scaled'}));
% Density is kg/m^3 and radii are km. Convert order-one gravity to SI:
Vz_SI = gravity{4} * 6.67430e-11 * 1000;
disp(table(Vz_SI,'VariableNames',{'Vz_m_per_s2'}));
magneticModel = [model(:,1:6) 0.5 0.2 -0.1];
magnetic = cell(1,10);
[magnetic{:}] = Tesseroid_MagneticEstimate(magneticModel,points,1e-12,1e-9);
assert(all(isfinite([gravity{:}]),'all') && all(isfinite([magnetic{:}]),'all'));
