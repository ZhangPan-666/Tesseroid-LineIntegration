% RunTimeRecord.m (TLI-MMSI port)
%
% Measured wall-clock times from the production TLI-MMSI cluster runs:
%
%   RunTime = [ grid_deg  time_all_20  time_Vz ];
%
% where
%   time_all_20  : wall-clock seconds for Work1_MassShellIntegrationMPI
%                  (all 20 gravitational quantities)
%   time_Vz      : wall-clock seconds for Work2_MassShellIntegrationVzMPI
%                  (Vz only)
%
% TLI-MMSI is a different physical discretization (grid-line integration):
% its absolute times and scaling are expected to differ from the 2-D MMSI
% values, so do NOT copy the original numbers.

% clc;clear;close all;

RunTime=...
[   10      1376.0     35.8;...
     5      1393.2     36.0;...
     2      1450.6     36.7;...
     1      1410.8     37.9;...
     0.5    1466.5     41.9;...
     0.25   1687.4     60.4;...
     0.1    1913.0     77.9;...
     0.05   2764.7     149.7];

% GridStep    20 Components    Vz
%
% These measured values are the values reported in the manuscript.
