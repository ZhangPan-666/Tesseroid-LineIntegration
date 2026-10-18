%% Plot Vz computed from ETOPO1 models with different grid spacings
% TLI-MMSI port of the original
%
% The ResultG.dat files are produced by the TLI-MMSI cluster drivers
%   TLI_MMSI/Work1_MassShellIntegrationMPI (all 20 quantities)
%   TLI_MMSI/Work2_MassShellIntegrationVzMPI (Vz only)
%
% The Work2 driver keeps the Work1 mass-result binary layout with
% V = Vx = Vy = NaN and Vz holding the computed values, so the Vz byte offset
% below (4 + 6*NumCoordinates*8) is identical to the original benchmark.

clc;
clear;
close all;

%% Input

GridSpacing=[10,5,2,1,0.5,0.25,0.1,0.05];
PanelLabel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

run('RunTimeRecord.m');
assert(exist('RunTime','var')==1 && size(RunTime,2)>=3, ...
    'RunTimeRecord.m must define RunTime as [grid, time20, timeVz].');
assert(all(isfinite(RunTime(:,2))) && all(isfinite(RunTime(:,3))), ...
    'RunTimeRecord.m must contain finite measured TLI-MMSI wall-clock times.');

%% Read the Vz component

PlotData=cell(size(GridSpacing));
LongitudeGrid=[];
LatitudeGrid=[];

for GridIndex=1:numel(GridSpacing)
    GridText=num2str(GridSpacing(GridIndex),'%g');
    FileName=['Results20G_' GridText 'Degree.dat'];

    [Longitude,Latitude,Vz]=ReadVzResult(FileName);

    if isempty(LongitudeGrid)
        LongitudeAxis=unique(Longitude,'stable');
        LatitudeAxis=unique(Latitude,'stable');
        NumLongitude=numel(LongitudeAxis);
        NumLatitude=numel(LatitudeAxis);

        assert(NumLongitude*NumLatitude==numel(Vz), ...
            'The observation coordinates do not form a regular grid.');

        LongitudeGrid=reshape(Longitude,NumLatitude,NumLongitude);
        LatitudeGrid=reshape(Latitude,NumLatitude,NumLongitude);
    else
        assert(numel(Vz)==numel(LongitudeGrid), ...
            'The result files use different observation grids.');
    end

    PlotData{GridIndex}=reshape(Vz,size(LongitudeGrid));
end

%% Shared colour range

MinimumValue=+inf;
MaximumValue=-inf;
for GridIndex=1:numel(PlotData)
    MinimumValue=min(MinimumValue, ...
        min(PlotData{GridIndex},[],'all','omitnan'));
    MaximumValue=max(MaximumValue, ...
        max(PlotData{GridIndex},[],'all','omitnan'));
end
ColourSpan=MaximumValue-MinimumValue;
if ~(isfinite(ColourSpan) && ColourSpan>0)
    ColourSpan=max(abs(MinimumValue),1);
end
ColourRange=[MinimumValue-0.02*ColourSpan, ...
    MaximumValue+0.02*ColourSpan];

ColourMap=MakeDivergingColourMap(256);

%% Work1/SAGP_Contour settings

GeoContourDirectory = getenv('TLI_GEOCONTOUR_DIR');
assert(isfolder(GeoContourDirectory), 'Set TLI_GEOCONTOUR_DIR to the external map plotting assets.');
addpath(GeoContourDirectory);

Coordinate.LonStart=min(LongitudeGrid,[],'all');
Coordinate.LonEnd=max(LongitudeGrid,[],'all');
Coordinate.LatStart=min(LatitudeGrid,[],'all');
Coordinate.LatEnd=max(LatitudeGrid,[],'all');
Coordinate.dLon=0.05;
Coordinate.dLat=0.05;

Contour.ContourColorMap=ColourMap;
Contour.ContourRange=ColourRange;

BoundaryDirectory=fullfile(GeoContourDirectory, ...
    'Coast and National Boundary');
Coast.CoastColor=[0,0,0];
Coast.CoastLineWidth=1.2;
Coast.CoastLineStyle='-';
Coast.CoastMain=shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L1.shp'));
Coast.LakesLine=shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L2.shp'));
Coast.CoastAntarctica=shaperead(fullfile(BoundaryDirectory, ...
    'GSHHS_l_L6.shp'));
Coast.NationalLine=shaperead(fullfile(BoundaryDirectory, ...
    'WDBII_border_f_L1.shp'));
TenDashData=load(fullfile(GeoContourDirectory, ...
    'Ten Dash Lines of China','TenDashLines.mat'),'TenDashLines');
Coast.TenDashLines=TenDashData.TenDashLines;

LonLatLine.LonLineInterval=5;
LonLatLine.LatLineInterval=5;
LonLatLine.LonLineColor=[0,0,0];
LonLatLine.LatLineColor=[0,0,0];
LonLatLine.LonLineWidth=0.8;
LonLatLine.LatLineWidth=0.8;
LonLatLine.LonLineStyle=':';
LonLatLine.LatLineStyle=':';
LonLatLine.TextShow=false;
LonLatLine.TextFont='Times New Roman';
LonLatLine.TextWeight='normal';
LonLatLine.TextFontSize=18;

%% Plot

F=figure('Color','white');
F.Units='pixels';
F.Position=[50,30,1650,1080];

% Runtime summary above the eight spatial panels.
RuntimeAxis=axes(F,'Position',[0.075,0.735,0.765,0.215]);
RuntimeAxis.FontName='Times New Roman';
RuntimeAxis.FontSize=18;
RuntimeAxis.LineWidth=1.0;
RuntimeAxis.XScale='log';
RuntimeAxis.XDir='reverse';
RuntimeAxis.XLim=[0.04,12];
RuntimeAxis.XTick=sort(GridSpacing);
RuntimeAxis.XTickLabel=compose('%g',sort(GridSpacing));
RuntimeAxis.XGrid='on';
RuntimeAxis.YGrid='on';
RuntimeAxis.GridAlpha=0.22;
RuntimeAxis.MinorGridAlpha=0.12;
hold(RuntimeAxis,'on');

yyaxis(RuntimeAxis,'left');
Runtime20Line=semilogx(RuntimeAxis,RunTime(:,1),RunTime(:,2),'-o', ...
    'Color',[0.10,0.36,0.68],'LineWidth',2.2,'MarkerSize',7, ...
    'MarkerFaceColor',[0.10,0.36,0.68]);
ylabel(RuntimeAxis,'Runtime for 20 quantities [s]', ...
    'FontName','Times New Roman','FontSize',20,'Color','black');
RuntimeAxis.YColor=[0.10,0.36,0.68];

yyaxis(RuntimeAxis,'right');
RuntimeVzLine=semilogx(RuntimeAxis,RunTime(:,1),RunTime(:,3),'-s', ...
    'Color',[0.78,0.20,0.16],'LineWidth',2.2,'MarkerSize',7, ...
    'MarkerFaceColor',[0.78,0.20,0.16]);
ylabel(RuntimeAxis,'Runtime for {\it V}_z [s]', ...
    'Interpreter','tex','FontName','Times New Roman', ...
    'FontSize',20,'Color','black');
RuntimeAxis.YColor=[0.78,0.20,0.16];

xlabel(RuntimeAxis,'Source-grid spacing [degree]', ...
    'FontName','Times New Roman','FontSize',20,'Color','black');
legend(RuntimeAxis,[Runtime20Line,RuntimeVzLine], ...
    {'20 quantities','{\it V}_z'},'Interpreter','tex', ...
    'FontName','Times New Roman','FontSize',17, ...
    'Location','northwest','Orientation','horizontal', ...
    'Box','off');

AxisWidth=0.185;
AxisHeight=0.285;
AxisLeft=[0.025,0.225,0.425,0.625];
AxisBottom=[0.390,0.055];
AxisHandles=gobjects(1,numel(GridSpacing));
AxisPositions=zeros(numel(GridSpacing),4);

for GridIndex=1:numel(GridSpacing)
    GridRow=ceil(GridIndex/4);
    GridColumn=mod(GridIndex-1,4)+1;
    AxisPositions(GridIndex,:)=[AxisLeft(GridColumn), ...
        AxisBottom(GridRow),AxisWidth,AxisHeight];
    AxisP=axes(F,'Position',AxisPositions(GridIndex,:));
    AxisP.PositionConstraint='innerposition';
    AxisHandles(GridIndex)=AxisP;

    Data.Longitude=LongitudeGrid(:);
    Data.Latitude=LatitudeGrid(:);
    Data.PlotData=PlotData{GridIndex}(:);
    LonLatLine.TextShow=(GridIndex==1);
    SAGP_Contour(AxisP,Data,Coordinate,Contour,Coast,LonLatLine);

    AxisP.XColor='none';
    AxisP.YColor='none';
    AxisP.FontName='Times New Roman';
    AxisP.FontSize=18;
    pbaspect(AxisP,[1,1,1]);

    text(AxisP,0.16,0.92,PanelLabel{GridIndex}, ...
        'Units','normalized','Interpreter','tex', ...
        'FontName','Times New Roman','FontSize',19, ...
        'FontWeight','normal','HorizontalAlignment','left', ...
        'VerticalAlignment','top','Color','black', ...
        'BackgroundColor','white','Margin',2);

end

ColourBar=colorbar(AxisHandles(end));
ColourBar.Position=[0.850,0.055,0.014,0.620];
for GridIndex=1:numel(AxisHandles)
    AxisHandles(GridIndex).Position=AxisPositions(GridIndex,:);
end
ColourBar.FontName='Times New Roman';
ColourBar.FontSize=18;
ColourBar.Label.String='Scaled {\it V}_z [m\cdot s^{-2}]';
ColourBar.Label.Interpreter='tex';
ColourBar.Label.FontName='Times New Roman';
ColourBar.Label.FontSize=22;

%% Export

drawnow;
exportgraphics(F,'RunTimeTest2_Vz.png', ...
    'Resolution',200,'BackgroundColor','white');
exportgraphics(F,'RunTimeTest2_Vz.jpg', ...
    'Resolution',200,'BackgroundColor','white');

%% Local functions

function [Longitude,Latitude,Vz]=ReadVzResult(FileName)
fid=fopen(FileName,'r','ieee-le');
assert(fid~=-1,'Failed to open result file: %s',FileName);
Cleaner=onCleanup(@() fclose(fid));

NumCoordinates=fread(fid,1,'int32');
assert(~isempty(NumCoordinates) && NumCoordinates>0, ...
    'Invalid coordinate count in %s.',FileName);
NumCoordinates=double(NumCoordinates);

Longitude=fread(fid,NumCoordinates,'double');
Latitude=fread(fid,NumCoordinates,'double');
assert(numel(Longitude)==NumCoordinates && ...
    numel(Latitude)==NumCoordinates, ...
    'Coordinate arrays are incomplete in %s.',FileName);

% Vz follows Longitude, Latitude, Radius, V, Vx and Vy.
% TLI Work2 writes V=Vx=Vy=NaN and Vz in the same Work1 layout.
VzOffset=4+6*NumCoordinates*8;
Status=fseek(fid,VzOffset,'bof');
assert(Status==0,'Failed to seek to Vz in %s.',FileName);
Vz=fread(fid,NumCoordinates,'double');
assert(numel(Vz)==NumCoordinates, ...
    'The Vz array is incomplete in %s.',FileName);
end

function ColourMap=MakeDivergingColourMap(NumberOfColours)
Position=[-100;-25;0;25;100];
Red=[0;135;255;240;230];
Green=[60;190;255;180;0];
Blue=[115;200;220;195;20];

QueryPosition=linspace(min(Position),max(Position),NumberOfColours);
ColourMap=[ ...
    interp1(Position,Red,QueryPosition,'linear')', ...
    interp1(Position,Green,QueryPosition,'linear')', ...
    interp1(Position,Blue,QueryPosition,'linear')']./255;
end
