%% write_MSIG_ForPar.m (TLI-MMSI port)
% Write the mass-shell parameter file required by
% Work2_MassShellIntegrationVzMPI.cpp in the TLI_MMSI directory.
%
% Binary layout (little-endian):
%   double  AbsTol, RelTol
%   double  Lons, Lone, Lats, Late, dLon, dLat
%   int32   NumLon, NumLat
%   double  R2(:), R1(:), rho(:)
%   int32   NumCoordinates
%   double  Longitude(:), Latitude(:), Radius(:)
%
% Convention:
%   R2, R1, rho are NumLat-by-NumLon matrices.
%   MATLAB column-major storage then matches the C++ indexing
%   idx = (k-1)*NumLat + l - 1, i.e. latitude varies fastest.
%
% This script is the TLI port of the original
% and is run once per source-grid spacing before the cluster run.
%
% The cluster driver is
% which reads MSIG.ForPar and writes ResultG.dat with the Work1 mass-result
% binary layout (V = Vx = Vy = NaN, Vz holds the computed values).
%
% Edit the USER INPUT section and run this script.

clear; clc;

SourceDir = getenv('TLI_ETOPO_DIR');
assert(isfolder(SourceDir), 'Set TLI_ETOPO_DIR to the directory containing ETOPO_<spacing>Degree.dat.');

%% ===================== USER INPUT =====================
outFile = 'MSIG.ForPar';
AbsTol  = 1e-16;
RelTol  = 1e-16;

% Grid centers
Lons = -180;
Lone =  180;
Lats = -90;
Late =  90;
dGrid = 0.05;

Lons=Lons+dGrid/2;
Lone=Lone-dGrid/2;
Lats=Lats+dGrid/2;
Late=Late-dGrid/2;

NumLon = round((Lone - Lons) / dGrid) + 1;
NumLat = round((Late - Lats) / dGrid) + 1;

dLon=dGrid;
dLat=dGrid;

% Model arrays: size must be [NumLat, NumLon]

Model=load(fullfile(SourceDir, ...
    ['ETOPO_' num2str(dGrid) 'Degree.dat']));
Model(:,3)=Model(:,3).*1e3;
a = 6378137;
f = 1.0 / 298.257223563; % WGS84 flattening
b=a-f*a;
e=sqrt(a^2-b^2)/a;

N=a./(sqrt(1-e^2.*sind(Model(:,2)).^2));

X=(N+Model(:,3)).*cosd(Model(:,2)).*cosd(Model(:,1));
Y=(N+Model(:,3)).*cosd(Model(:,2)).*sind(Model(:,1));
Z=(N*(1-e^2)+Model(:,3)).*sind(Model(:,2));

[~,~,R2]=cart2sph(X,Y,Z);
R2=R2./1e3;

R1Altitude=-20000;
X=(N+R1Altitude).*cosd(Model(:,2)).*cosd(Model(:,1));
Y=(N+R1Altitude).*cosd(Model(:,2)).*sind(Model(:,1));
Z=(N*(1-e^2)+R1Altitude).*sind(Model(:,2));

[~,~,R1]=cart2sph(X,Y,Z);
R1=R1./1e3;

R2=reshape(R2,NumLat,NumLon);
R1=reshape(R1,NumLat,NumLon);
rho = 1 * ones(NumLat, NumLon);

% Computational points (can be row or column vectors)
[Longitude,Latitude]=meshgrid(130:0.05:150,10:0.05:40);
Longitude=Longitude(:);
Latitude=Latitude(:);
Radius=(6371.2+15).*ones(size(Longitude));

%% ======================================================

% Checks
assert(all(size(R1)  == [NumLat, NumLon]), 'R1 must be NumLat-by-NumLon.');
assert(all(size(R2)  == [NumLat, NumLon]), 'R2 must be NumLat-by-NumLon.');
assert(all(size(rho) == [NumLat, NumLon]), 'rho must be NumLat-by-NumLon.');

Longitude = Longitude(:);
Latitude  = Latitude(:);
Radius    = Radius(:);
NumCoordinates = numel(Longitude);
assert(numel(Latitude) == NumCoordinates && numel(Radius) == NumCoordinates, ...
    'Longitude, Latitude, Radius must have the same length.');

fid = fopen(outFile, 'w', 'ieee-le');
assert(fid ~= -1, 'Failed to open output file: %s', outFile);
cleaner = onCleanup(@() fclose(fid));

fwrite(fid, double(AbsTol), 'double');
fwrite(fid, double(RelTol), 'double');

fwrite(fid, double(Lons), 'double');
fwrite(fid, double(Lone), 'double');
fwrite(fid, double(Lats), 'double');
fwrite(fid, double(Late), 'double');
fwrite(fid, double(dLon), 'double');
fwrite(fid, double(dLat), 'double');

fwrite(fid, int32(NumLon), 'int32');
fwrite(fid, int32(NumLat), 'int32');

fwrite(fid, double(R2(:)),  'double');
fwrite(fid, double(R1(:)),  'double');
fwrite(fid, double(rho(:)), 'double');

fwrite(fid, int32(NumCoordinates), 'int32');
fwrite(fid, double(Longitude), 'double');
fwrite(fid, double(Latitude),  'double');
fwrite(fid, double(Radius),    'double');

disp(['Mass-shell parameter file written: ', outFile]);
