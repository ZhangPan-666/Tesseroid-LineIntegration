%%

clc;
clear;
close all;

%%

FileName='../RuntimeTest_LaplacePoisson_HeightPlus1km.csv';

Data=readmatrix(FileName);

Longitude=Data(:,2);
Latitude=Data(:,3);
Runtime=Data(:,8:27);

LongitudeAxis=unique(Longitude,'stable');
LatitudeAxis=unique(Latitude,'stable');

NumLongitude=length(LongitudeAxis);
NumLatitude=length(LatitudeAxis);

if NumLongitude*NumLatitude~=size(Data,1)
    error('The CSV file does not contain a complete longitude-latitude grid.');
end

RuntimeGrid=zeros(NumLatitude,NumLongitude,20);
for ComponentIndex=1:20
    RuntimeGrid(:,:,ComponentIndex)=reshape( ...
        Runtime(:,ComponentIndex),NumLongitude,NumLatitude).';
end

%%

ComponentNames={ ...
    '\it{V}', ...
    '\it{Vx}','\it{Vy}','\it{Vz}', ...
    '\it{Vxx}','\it{Vxy}','\it{Vyy}','\it{Vzx}', ...
    '\it{Vzy}','\it{Vzz}', ...
    '\it{Vxxx}','\it{Vxxy}','\it{Vxxz}','\it{Vxyz}', ...
    '\it{Vyyx}','\it{Vyyy}','\it{Vyyz}','\it{Vzzx}', ...
    '\it{Vzzy}','\it{Vzzz}'};

PanelNames=arrayfun(@(Index)char('a'+Index-1),1:20, ...
    'UniformOutput',false);

PositiveRuntime=Runtime(Runtime>0&isfinite(Runtime));
if isempty(PositiveRuntime)
    error('No positive finite computation times were found in the CSV file.');
end

ColorLimits=10.^[floor(log10(min(PositiveRuntime))), ...
    ceil(log10(max(PositiveRuntime)))];

if ColorLimits(1)==ColorLimits(2)
    ColorLimits=ColorLimits.*[0.1 10];
end

ColorTicks=10.^(log10(ColorLimits(1)):log10(ColorLimits(2)));

%%

ColorAnchor=[ ...
      0   0   4; ...
     22  11  57; ...
     66  10 104; ...
    106  23 110; ...
    147  38 103; ...
    188  55  84; ...
    221  81  58; ...
    243 120  25; ...
    252 170  53; ...
    252 220 120; ...
    252 255 164]./255;

RuntimeColorMap=interp1(linspace(0,1,size(ColorAnchor,1)), ...
    ColorAnchor,linspace(0,1,256),'linear');

%%

F=figure('Visible','off');
F.Units='pixels';
F.Position=[50,50,1450,1800];
F.Color=[1 1 1];
F.InvertHardcopy='off';
F.PaperUnits='inches';
F.PaperPosition=[0 0 14.5 18];
F.PaperSize=[14.5 18];
F.PaperPositionMode='manual';

LeftPosition=0.055;
BottomPosition=0.055;
AxisWidth=0.195;
AxisHeight=0.157;
HorizontalGap=0.010;
VerticalGap=0.025;

XTick=[105 120 135];
YTick=[25 40 55];

for ComponentIndex=1:20
    RowIndex=ceil(ComponentIndex/4);
    ColumnIndex=mod(ComponentIndex-1,4)+1;

    AxisP=axes;
    AxisP.Position=[ ...
        LeftPosition+(ColumnIndex-1)*(AxisWidth+HorizontalGap), ...
        BottomPosition+(5-RowIndex)*(AxisHeight+VerticalGap), ...
        AxisWidth,AxisHeight];

    imagesc(AxisP,LongitudeAxis,LatitudeAxis, ...
        RuntimeGrid(:,:,ComponentIndex));

    AxisP.FontName='times new roman';
    AxisP.FontSize=22;
    AxisP.FontSizeMode='manual';
    AxisP.LineWidth=1.5;
    AxisP.Box='on';

    AxisP.XTick=XTick;
    AxisP.YTick=YTick;
    AxisP.TickDir='out';
    AxisP.TickDirMode='manual';
    AxisP.YDir='normal';
    AxisP.ColorScale='log';
    AxisP.CLim=ColorLimits;
    axis(AxisP,'equal');
    AxisP.XLim=[min(LongitudeAxis) max(LongitudeAxis)];
    AxisP.YLim=[min(LatitudeAxis) max(LatitudeAxis)];

    if ColumnIndex==1
        AxisP.YLabel.String=['Latitude (' char(176) 'N)'];
        AxisP.YLabel.Interpreter='none';
    else
        AxisP.YTickLabel=[];
    end

    if RowIndex==5
        AxisP.XLabel.String=['Longitude (' char(176) 'E)'];
        AxisP.XLabel.Interpreter='none';
    else
        AxisP.XTickLabel=[];
    end

    TP=title(AxisP,['(' PanelNames{ComponentIndex} ') ' ...
        ComponentNames{ComponentIndex}]);
    TP.FontName='times new roman';
    TP.FontWeight='normal';
    TP.FontSize=23;
    TP.Color=[0 0 0];
    TP.Interpreter='tex';
    TP.Units='normalized';
    TP.Position=[0 1.015 0];
    TP.HorizontalAlignment='left';

    ComponentRuntime=Runtime(:,ComponentIndex);
    ComponentRuntime=ComponentRuntime( ...
        ComponentRuntime>0&isfinite(ComponentRuntime));

    MinRuntime=min(ComponentRuntime);
    MaxRuntime=max(ComponentRuntime);

    MinExponent=floor(log10(MinRuntime));
    MaxExponent=floor(log10(MaxRuntime));
    MinMantissa=MinRuntime/10^MinExponent;
    MaxMantissa=MaxRuntime/10^MaxExponent;

    MinMantissa=round(MinMantissa,2);
    MaxMantissa=round(MaxMantissa,2);

    if MinMantissa>=10
        MinMantissa=MinMantissa/10;
        MinExponent=MinExponent+1;
    end

    if MaxMantissa>=10
        MaxMantissa=MaxMantissa/10;
        MaxExponent=MaxExponent+1;
    end

    SP=text(AxisP,0.025,0.965,sprintf( ...
        ['min = %.2f\\times10^{%d} s\n' ...
        'max = %.2f\\times10^{%d} s'], ...
        MinMantissa,MinExponent,MaxMantissa,MaxExponent));
    SP.Units='normalized';
    SP.Position=[0.025 0.965 0];
    SP.FontName='times new roman';
    SP.FontWeight='normal';
    SP.FontSize=18;
    SP.Color=[1 1 1];
    SP.Interpreter='tex';
    SP.HorizontalAlignment='left';
    SP.VerticalAlignment='top';

    colormap(AxisP,RuntimeColorMap);
end

%%

c=colorbar(AxisP);
c.Position=[0.878 BottomPosition 0.018 ...
    5*AxisHeight+4*VerticalGap];
c.Limits=ColorLimits;
c.Ticks=ColorTicks;
c.FontName='times new roman';
c.FontSize=22;
c.Label.String='Computation time [s]';
c.Label.FontName='times new roman';
c.Label.FontSize=24;
c.Label.VerticalAlignment='middle';

%%

drawnow;
exportgraphics(F,'RuntimeTest_1km.png', ...
    'Resolution',600,'BackgroundColor','white');
RuntimeImage=imread('RuntimeTest_1km.png');
imwrite(RuntimeImage,'RuntimeTest_1km.jpg','Quality',95);
