function build_Tesseroid_LineIntegration_mex(useOpenMP)
%BUILD_TESSEROID_LINEINTEGRATION_MEX Build with MATLAB's selected compiler.
%   Requires mex -setup C++. Optional logical argument enables OpenMP.
%   Default: true on Windows/Linux, false on macOS.
if nargin < 1, useOpenMP = ~ismac; end
validateattributes(useOpenMP, {'logical','numeric'}, {'scalar'});
functionsDirectory = fileparts(mfilename('fullpath'));
includeDirectory = fullfile(fileparts(functionsDirectory),'C++Codes');
source = fullfile(functionsDirectory,'tli_lineintegration_mex.cpp');
args = {'-R2018a','-O',['-I' includeDirectory],'-outdir',functionsDirectory,source};
if ispc
    flags = 'COMPFLAGS=$COMPFLAGS /std:c++17 /EHsc';
    if useOpenMP, flags = [flags ' /openmp']; end
    args = [{flags},args];
else
    flags = 'CXXFLAGS=$CXXFLAGS -std=c++17';
    if useOpenMP
        flags = [flags ' -fopenmp'];
        args = [{'LDFLAGS=$LDFLAGS -fopenmp'},args];
    end
    args = [{flags},args];
end
clear tli_lineintegration_mex
mex(args{:});
rehash
fprintf('Built %s (OpenMP requested: %d)\n', ...
    fullfile(functionsDirectory,['tli_lineintegration_mex.' mexext]),useOpenMP);
end
