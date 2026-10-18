%% write_MSIM_ForPar.m (TLI-MMSI port)
% Write the magnetization-shell parameter file required by
% Work1_MagnetizationShellIntegrationMPI.cpp (all 10 magnetic quantities)
% and Work2_MagnetizationShellIntegrationVectorMPI.cpp (Vx, Vy, Vz only)
% in the TLI_MMSI directory.
%
% Binary layout (little-endian):
%   double  AbsTol, RelTol
%   double  Lons, Lone, Lats, Late, dLon, dLat
%   int32   NumLon, NumLat
%   double  R2(:), R1(:), Mx(:), My(:), Mz(:)
%   int32   NumCoordinates
%   double  Longitude(:), Latitude(:), Radius(:)
%
% Convention:
%   R2, R1, Mx, My, Mz are NumLat-by-NumLon matrices.
%   MATLAB column-major storage then matches the C++ indexing
%   idx = (k-1)*NumLat + l - 1, i.e. latitude varies fastest.
%
% This script is the TLI-MMSI port of the original
%
% Edit the USER INPUT section and run this script.

clear; clc;
here = fileparts(mfilename('fullpath'));
vimFile = getenv('TLI_VIM_FILE');
assert(isfile(vimFile), 'Set TLI_VIM_FILE to the external 7-VIMModel.dat file.');
VIMModel=readmatrix(vimFile);

%% ===================== USER INPUT =====================
outFile = 'MSIM.ForPar';
AbsTol  = 1e-16;
RelTol  = 1e-16;

% Grid centers
Lons = -180.0;
Lone =  180.0;
Lats = -90.0;
Late =  90.0;
dLon =  0.25;
dLat =  0.25;

Lons=Lons+dLon/2;
Lone=Lone-dLon/2;
Lats=Lats+dLat/2;
Late=Late-dLat/2;

NumLon = round((Lone - Lons) / dLon) + 1;
NumLat = round((Late - Lats) / dLat) + 1;

% Model arrays: size must be [NumLat, NumLon]
R1 = reshape(VIMModel(:,5),NumLat,NumLon);
R2 = reshape(VIMModel(:,6),NumLat,NumLon);
Mx = reshape(VIMModel(:,7),NumLat,NumLon);
My = reshape(VIMModel(:,8),NumLat,NumLon);
Mz = reshape(VIMModel(:,9),NumLat,NumLon);

% Computational points (can be row or column vectors)
[Longitude,Latitude]=meshgrid(130:0.05:150,10:0.05:40);
Longitude=Longitude(:);
Latitude=Latitude(:);
Radius=(6371.2+15).*ones(size(Longitude));

%% ======================================================

% Checks
assert(all(size(R1) == [NumLat, NumLon]), 'R1 must be NumLat-by-NumLon.');
assert(all(size(R2) == [NumLat, NumLon]), 'R2 must be NumLat-by-NumLon.');
assert(all(size(Mx) == [NumLat, NumLon]), 'Mx must be NumLat-by-NumLon.');
assert(all(size(My) == [NumLat, NumLon]), 'My must be NumLat-by-NumLon.');
assert(all(size(Mz) == [NumLat, NumLon]), 'Mz must be NumLat-by-NumLon.');

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

fwrite(fid, double(R2(:)), 'double');
fwrite(fid, double(R1(:)), 'double');
fwrite(fid, double(Mx(:)), 'double');
fwrite(fid, double(My(:)), 'double');
fwrite(fid, double(Mz(:)), 'double');

fwrite(fid, int32(NumCoordinates), 'int32');
fwrite(fid, double(Longitude), 'double');
fwrite(fid, double(Latitude),  'double');
fwrite(fid, double(Radius),    'double');

disp(['Magnetization-shell parameter file written: ', outFile]);
