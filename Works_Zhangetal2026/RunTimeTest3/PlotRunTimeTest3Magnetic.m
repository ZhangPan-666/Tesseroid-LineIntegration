%% Plot magnetic potential, field, and gradient components
% TLI-MMSI port of the original
%
% The ResultM.dat file is produced by the TLI-MMSI cluster drivers
%   TLI_MMSI/Work1_MagnetizationShellIntegrationMPI (all 10 magnetic
%   quantities)
%   TLI_MMSI/Work2_MagnetizationShellIntegrationVectorMPI (Vx, Vy, Vz only)
%
% The Work2 driver keeps the Work1 magnetic-result binary layout
% (Longitude, Latitude, Radius, V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy, Vzz)
% with V = Vxx = Vxy = Vyy = Vzx = Vzy = Vzz = 0.

clc;
clear;
close all;

%% Input

InputFile = 'ResultM.dat';
OutputName = 'fig9';

ComponentNames = {'V','Vx','Vy','Vz','Vxx','Vxy','Vyy','Vzx','Vzy','Vzz'};
PanelLabels = {'(a)','(b)','(c)','(d)','(e)', ...
    '(f)','(g)','(h)','(i)','(j)'};
UnitLabels = {'mT\cdot m','nT','nT','nT', ...
    'pT/m','pT/m','pT/m','pT/m','pT/m','pT/m'};
ConversionFactors = [0.1, repmat(100, 1, 9)];

run('RunTimeRecord.m');
assert(exist('RunTime','var')==1 && size(RunTime,2)>=3, ...
    'RunTimeRecord.m must define RunTime as [grid, timeAll10, timeVector].');
assert(all(isfinite(RunTime(:,2))) && all(isfinite(RunTime(:,3))), ...
    'RunTimeRecord.m must contain finite measured TLI-MMSI wall-clock times.');

M = readMSIMResults(InputFile);

LongitudeAxis = unique(M.Longitude, 'stable');
LatitudeAxis = unique(M.Latitude, 'stable');
numLongitude = numel(LongitudeAxis);
numLatitude = numel(LatitudeAxis);
assert(numLongitude*numLatitude == M.NumCoordinates, ...
    'The observation coordinates do not form a regular grid.');

LongitudeGrid = reshape(M.Longitude, numLatitude, numLongitude);
LatitudeGrid = reshape(M.Latitude, numLatitude, numLongitude);
PlotData = cell(1, numel(ComponentNames));
ColourRanges = zeros(numel(ComponentNames), 2);

for componentIndex = 1:numel(ComponentNames)
    % C++ returns the kernel combinations before mu0/(4*pi), with radii in km.
    values = M.(ComponentNames{componentIndex}) .* ...
        ConversionFactors(componentIndex);
    PlotData{componentIndex} = reshape(values, numLatitude, numLongitude);
    ColourRanges(componentIndex, :) = componentColourRange(values);
end

% Use a fixed potential range that retains the full computed extrema while
% improving contrast across the predominantly negative regional field.
ColourRanges(1, :) = [-8, 1];

ColourMap = makeDivergingColourMap(256);

%% Work1/SAGP_Contour settings

GeoContourDirectory = getenv('TLI_GEOCONTOUR_DIR');
assert(isfolder(GeoContourDirectory), 'Set TLI_GEOCONTOUR_DIR to the external map plotting assets.');
addpath(GeoContourDirectory);

Coordinate.LonStart = min(LongitudeGrid, [], 'all');
Coordinate.LonEnd = max(LongitudeGrid, [], 'all');
Coordinate.LatStart = min(LatitudeGrid, [], 'all');
Coordinate.LatEnd = max(LatitudeGrid, [], 'all');
Coordinate.dLon = median(diff(LongitudeAxis));
Coordinate.dLat = median(diff(LatitudeAxis));

Contour.ContourColorMap = ColourMap;

BoundaryDirectory = fullfile(GeoContourDirectory, ...
    'Coast and National Boundary');
Coast.CoastColor = [0,0,0];
Coast.CoastLineWidth = 1.1;
Coast.CoastLineStyle = '-';
Coast.CoastMain = shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L1.shp'));
Coast.LakesLine = shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L2.shp'));
Coast.CoastAntarctica = shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L6.shp'));
Coast.NationalLine = shaperead(fullfile(BoundaryDirectory, ...
    'WDBII_border_f_L1.shp'));
tenDashData = load(fullfile(GeoContourDirectory, ...
    'Ten Dash Lines of China','TenDashLines.mat'), 'TenDashLines');
Coast.TenDashLines = tenDashData.TenDashLines;

LonLatLine.LonLineInterval = 5;
LonLatLine.LatLineInterval = 5;
LonLatLine.LonLineColor = [0,0,0];
LonLatLine.LatLineColor = [0,0,0];
LonLatLine.LonLineWidth = 0.7;
LonLatLine.LatLineWidth = 0.7;
LonLatLine.LonLineStyle = ':';
LonLatLine.LatLineStyle = ':';
LonLatLine.TextShow = false;
LonLatLine.TextFont = 'Times New Roman';
LonLatLine.TextWeight = 'normal';
LonLatLine.TextFontSize = 22;

%% Plot

F = figure('Color', 'white', 'Visible', 'off');
F.Units = 'pixels';
F.Position = [50,50,2300,1050];

axisWidth = 0.145;
axisHeight = 0.39;
axisLeft = [0.022,0.214,0.406,0.598,0.790];
axisBottom = [0.515,0.170];
plotBoxHeight = axisWidth * F.Position(3) / F.Position(4);
colourBarHeight = 0.84 * plotBoxHeight;
ColourBars = gobjects(1, numel(ComponentNames));

for componentIndex = 1:numel(ComponentNames)
    row = ceil(componentIndex/5);
    column = mod(componentIndex-1, 5) + 1;
    axisPosition = [axisLeft(column),axisBottom(row),axisWidth,axisHeight];
    AxisP = axes(F, 'Position', axisPosition);
    AxisP.PositionConstraint = 'innerposition';

    Data.Longitude = LongitudeGrid(:);
    Data.Latitude = LatitudeGrid(:);
    Data.PlotData = PlotData{componentIndex}(:);
    Contour.ContourRange = ColourRanges(componentIndex, :);
    LonLatLine.TextShow = componentIndex == 1;
    SAGP_Contour(AxisP, Data, Coordinate, Contour, Coast, LonLatLine);

    AxisP.XColor = 'none';
    AxisP.YColor = 'none';
    AxisP.FontName = 'Times New Roman';
    AxisP.FontSize = 22;
    pbaspect(AxisP, [1,1,1]);

    text(AxisP, 0.16, 0.92, PanelLabels{componentIndex}, ...
        'Units', 'normalized', 'Interpreter', 'tex', ...
        'FontName', 'Times New Roman', 'FontSize', 24, ...
        'FontWeight', 'normal', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Color', 'black', ...
        'BackgroundColor', 'white', 'Margin', 2);

    ColourBar = colorbar(AxisP);
    plotBoxBottom = axisBottom(row) + (axisHeight-plotBoxHeight)/2;
    colourBarBottom = plotBoxBottom + (plotBoxHeight-colourBarHeight)/2;
    ColourBar.Position = [axisLeft(column)+0.151, ...
        colourBarBottom,0.008,colourBarHeight];
    AxisP.Position = axisPosition;
    ColourBar.FontName = 'Times New Roman';
    ColourBar.FontSize = 15;
    ColourBar.Label.String = UnitLabels{componentIndex};
    ColourBar.Label.Interpreter = 'tex';
    ColourBar.Label.FontName = 'Times New Roman';
    ColourBar.Label.FontSize = 19;
    if componentIndex == 1
        ColourBar.Ticks = [-8,-6,-4,-2,0,1];
    end
    ColourBars(componentIndex) = ColourBar;
end

% Later axes must not obscure colour-bar ticks or labels from earlier panels.
uistack(ColourBars, 'top');

%% Export

drawnow;
exportgraphics(F, [OutputName '.png'], ...
    'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(F, [OutputName '.jpg'], ...
    'Resolution', 600, 'BackgroundColor', 'white');
close(F);

%% Local functions

function M = readMSIMResults(filename)
fid = fopen(filename, 'r', 'ieee-le');
assert(fid ~= -1, 'Failed to open result file: %s', filename);
cleaner = onCleanup(@() fclose(fid));

M = struct();
M.NumCoordinates = double(fread(fid, 1, 'int32'));
assert(~isempty(M.NumCoordinates) && M.NumCoordinates > 0, ...
    'Failed to read NumCoordinates from %s.', filename);

fieldNames = {'Longitude','Latitude','Radius','V','Vx','Vy','Vz', ...
    'Vxx','Vxy','Vyy','Vzx','Vzy','Vzz'};
for fieldIndex = 1:numel(fieldNames)
    M.(fieldNames{fieldIndex}) = fread(fid, M.NumCoordinates, 'double');
    assert(numel(M.(fieldNames{fieldIndex})) == M.NumCoordinates, ...
        'Incomplete %s array in %s.', fieldNames{fieldIndex}, filename);
end
end

function colourRange = componentColourRange(values)
minimumValue = min(values, [], 'all', 'omitnan');
maximumValue = max(values, [], 'all', 'omitnan');
assert(isfinite(minimumValue) && isfinite(maximumValue), ...
    'Component contains no finite values.');

if minimumValue < 0 && maximumValue > 0
    limit = max(abs([minimumValue, maximumValue]));
    colourRange = [-limit, limit];
else
    span = maximumValue - minimumValue;
    if span == 0
        span = max(abs(minimumValue), 1);
    end
    colourRange = [minimumValue-0.02*span, maximumValue+0.02*span];
end
end

function colourMap = makeDivergingColourMap(numberOfColours)
position = [-100;-25;0;25;100];
red = [0;135;255;240;230];
green = [60;190;255;180;0];
blue = [115;200;220;195;20];

queryPosition = linspace(min(position), max(position), numberOfColours);
colourMap = [ ...
    interp1(position, red, queryPosition, 'linear')', ...
    interp1(position, green, queryPosition, 'linear')', ...
    interp1(position, blue, queryPosition, 'linear')'] ./ 255;
end
