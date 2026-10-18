% Test No. 6: R_E - 31 km (inside the spherical shell, below it).
clc; clear; close all;
here = fileparts(mfilename('fullpath'));
dataFile = fullfile(here, '..', 'results', 'ResultsPan_10x10Split_-31km.dat');
outputBase = fullfile(here, 'Results2_TLImmSI_-31km');
plot_shell_figure(dataFile, outputBase);
fprintf('Done: %s.png / %s.jpg\n', outputBase, outputBase);
