% Test No. 7: R_E - 15 km (within the spherical shell).
% Reproduces the Figure4 'Pan -15km' panel style for the TLI-MMSI result.

clc; clear; close all;

here = fileparts(mfilename('fullpath'));
dataFile = fullfile(here, '..', 'results', 'ResultsPan_10x10Split_-15km.dat');
outputBase = fullfile(here, 'Results3_TLImmSI_-15km');

plot_shell_figure(dataFile, outputBase);

fprintf('Done: %s.png / %s.jpg\n', outputBase);