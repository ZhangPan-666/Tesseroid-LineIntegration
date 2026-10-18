% RunTimeRecord.m (TLI-MMSI port)
%
% Measured wall-clock times from the production TLI-MMSI magnetization runs:
%
%   RunTime = [ grid_deg  time_all_10  time_vector ];
%
% where
%   time_all_10  : wall-clock seconds for Work1_MagnetizationShellIntegrationMPI
%                  (all 10 magnetic quantities)
%   time_vector  : wall-clock seconds for Work2_MagnetizationShellIntegrationVectorMPI
%                  (the three-component magnetic-field vector Vx, Vy, Vz)
%
% The original Tesseroid RunTimeTest3 figures for the 2-D MMSI method were
%
%   grid 0.25: 34754.6 s (10 quantities) / 12958.8 s (vector)
%
% TLI-MMSI is a different physical discretization (grid-line integration):
% its absolute times and the all-10 / vector ratio are expected to differ
% from the 2-D MMSI values, so do NOT copy the original numbers.

% clc;clear;close all;

RunTime=...
[0.25    3419.4     971.5];

% GridStep    10 Components    Vector
%
% These measured values are the values reported in the manuscript.
