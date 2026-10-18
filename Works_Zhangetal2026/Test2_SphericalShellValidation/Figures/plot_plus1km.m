% Test No. 5: R_E + 1 km (outside the spherical shell).
clc; clear; close all;
here = fileparts(mfilename('fullpath'));
dataFile = fullfile(here, '..', 'results', 'ResultsPan_10x10Split_1km.dat');
outputBase = fullfile(here, 'Results1_TLImmSI_1km');
plot_shell_figure(dataFile, outputBase);
fprintf('Done: %s.png / %s.jpg\n', outputBase, outputBase);
