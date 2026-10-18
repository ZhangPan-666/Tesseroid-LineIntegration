%% Compatibility entry point for the complete runtime audit
% The authoritative analysis is maintained in Analyze_RuntimeTest1.m.

here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'Analyze_RuntimeTest1.m'));
