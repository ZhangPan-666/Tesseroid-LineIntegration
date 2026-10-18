%%

clc;clear;close all;

SourceDir = getenv('TLI_ETOPO_DIR');
assert(isfolder(SourceDir), 'Set TLI_ETOPO_DIR to the directory containing ETOPO_<spacing>Degree.dat.');

%%
% Generate the global cell-centred coordinate grid for a given grid spacing.
% This is the TLI-MMSI port of the original
% The output Coordinate_<dgrid>Degree.dat is used together with the
% ETOPO_<dgrid>Degree.dat files for RunTimeTest_Gravity.m.

Lons = -180;
Lone =  180;
Lats = -90;
Late =  90;
dgrid =  0.1;

Lons=Lons+dgrid/2;
Lone=Lone-dgrid/2;
Lats=Lats+dgrid/2;
Late=Late-dgrid/2;

%%

[Longitude,Latitude]=meshgrid(Lons:dgrid:Lone,Lats:dgrid:Late);
Longitude=Longitude(:);
Latitude=Latitude(:);

writematrix([Longitude,Latitude],fullfile(SourceDir, ...
    ['Coordinate_' num2str(dgrid) 'Degree.dat']))
