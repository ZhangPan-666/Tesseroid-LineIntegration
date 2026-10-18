function PlotTest4FourRows(Longitude,Latitude,PlotData, ...
    SecondLimits,ThirdXLimits,ThirdYLimits,ThirdZLimits,OutputName)

%%

if ~isequal(size(PlotData),[4,5])
    error('PlotData must be a 4-by-5 cell array.');
end

LongitudeAxis=unique(Longitude(:));
LatitudeAxis=unique(Latitude(:));

%%

cbar1=MakeColorBar([-100;-10;0;10;100],...Position
    [0;135;255;240;230],...Red
    [60;190;255;180;0],...Green
    [115;200;220;195;20]);% Blue

cbar2=MakeColorBar([-18;-15;-12;-9;-6;-3;0],...Position
    [45;55;77;152;255;255;255],...Red
    [45;126;175;78;127;255;255],...Green
    [45;184;74;163;0;51;255]);% Blue

%%

ComponentTitles={ ...
    '\it{Vxx}','\it{Vyy}','\it{Vzz}','\it{Vxx+Vyy+Vzz}'; ...
    '\it{Vxxx}','\it{Vyyx}','\it{Vzzx}','\it{Vxxx+Vyyx+Vzzx}'; ...
    '\it{Vxxy}','\it{Vyyy}','\it{Vzzy}','\it{Vxxy+Vyyy+Vzzy}'; ...
    '\it{Vxxz}','\it{Vyyz}','\it{Vzzz}','\it{Vxxz+Vyyz+Vzzz}'};

%%

F=figure('Visible','off');
F.Units='pixels';
F.Position=[50,50,1600,1300];
F.Color=[1 1 1];

LeftPosition=[0.045 0.197 0.349 0.501 0.728];
BottomPosition=0.09;
AxisWidth=0.145;
AxisHeight=0.1785;
VerticalGap=0.052;

XTick=[105 120 135];
YTick=[25 40 55];

PanelLetters={'a','b','c','d','e','f','g','h','i','j', ...
             'k','l','m','n','o','p','q','r','s','t'};

for RowIndex=1:4
    RowBottom=BottomPosition+(4-RowIndex)*(AxisHeight+VerticalGap);
    GroupIndex=RowIndex;

    if GroupIndex==1
        DataLimits=SecondLimits;
        UnitLabel='Scaled [s^{-2}]';
    elseif GroupIndex==2
        DataLimits=ThirdXLimits;
        UnitLabel='Scaled [m^{-1}\cdot s^{-2}]';
    elseif GroupIndex==3
        DataLimits=ThirdYLimits;
        UnitLabel='Scaled [m^{-1}\cdot s^{-2}]';
    else
        DataLimits=ThirdZLimits;
        UnitLabel='Scaled [m^{-1}\cdot s^{-2}]';
    end

    for ColumnIndex=1:5
        AxisP=axes;
        AxisP.Position=[LeftPosition(ColumnIndex),RowBottom, ...
            AxisWidth,AxisHeight];

        Value=PlotData{RowIndex,ColumnIndex};
        imagesc(AxisP,LongitudeAxis,LatitudeAxis,Value);

        if ColumnIndex<=4
            AxisP.Children.CData=GetCData(Value,cbar1,DataLimits);
            AxisP.CLim=DataLimits;
            colormap(AxisP,cbar1);
        else
            AxisP.Children.CData=GetCData(Value,cbar2,[-18,0]);
            AxisP.CLim=[-18,0];
            colormap(AxisP,cbar2);
        end

        AxisP.FontName='times new roman';
        AxisP.FontSize=19;
        AxisP.FontSizeMode='manual';
        AxisP.LineWidth=1.3;
        AxisP.Box='on';
        AxisP.XTick=XTick;
        AxisP.XTickLabelRotation=0;
        AxisP.YTick=YTick;
        AxisP.TickDir='out';
        AxisP.TickDirMode='manual';
        AxisP.YDir='normal';
        AxisP.XLim=[100,140];
        AxisP.YLim=[20,60];
        axis(AxisP,'equal');

        if ColumnIndex==1
            AxisP.YLabel.String='Latitude (^o)';
            AxisP.YLabel.FontSize=22;
        else
            AxisP.YTickLabel=[];
        end

        if RowIndex==4
            AxisP.XLabel.String='Longitude (^o)';
            AxisP.XLabel.FontSize=22;
        else
            AxisP.XTickLabel=[];
        end

        PanelLettersIndex=(RowIndex-1)*5+ColumnIndex;
        PanelName=PanelLetters{PanelLettersIndex};

        if ColumnIndex<=4
            TitleString=['(' PanelName ') ' ...
                ComponentTitles{GroupIndex,ColumnIndex}];
        else
            TitleString=['(' PanelName ') Absolute error'];
        end

        TP=title(AxisP,TitleString);
        TP.FontName='times new roman';
        TP.FontWeight='normal';
        TP.FontSize=20;
        TP.Interpreter='tex';
        TP.Units='normalized';
        TP.Position=[0 1.01 0];
        TP.HorizontalAlignment='left';
    end

    DataAxis=findobj(F,'Type','axes');
    DataAxis=DataAxis(2);

    c=colorbar(DataAxis);
    c.Position=[0.655 RowBottom 0.011 AxisHeight];
    c.Limits=DataLimits;
    c.Ticks=[DataLimits(1),0,DataLimits(2)];
    c.FontName='times new roman';
    c.FontSize=16;
    c.Label.String=UnitLabel;
    c.Label.FontName='times new roman';
    c.Label.FontSize=20;
    c.Label.Interpreter='tex';
    c.Label.Units='normalized';
    c.Label.Position=[3.6 0.5 0];
    c.Label.Rotation=90;
    c.Label.HorizontalAlignment='center';
    c.Label.VerticalAlignment='middle';

    ErrorAxis=findobj(F,'Type','axes');
    ErrorAxis=ErrorAxis(1);

    c=colorbar(ErrorAxis);
    c.Position=[0.882 RowBottom 0.011 AxisHeight];
    c.Limits=[-18,0];
    c.Ticks=[-18 -15 -12 -9 -6 -3 0];
    c.FontName='times new roman';
    c.FontSize=16;
    c.Label.String='log_{10} |scaled error|';
    c.Label.FontName='times new roman';
    c.Label.FontSize=20;
    c.Label.Interpreter='tex';
    c.Label.Units='normalized';
    c.Label.Position=[5.2 0.5 0];
    c.Label.Rotation=90;
    c.Label.HorizontalAlignment='center';
    c.Label.VerticalAlignment='middle';
end

%%

drawnow;
exportgraphics(F,[OutputName '.png'], ...
    'Resolution',600,'BackgroundColor','white');
OutputImage=imread([OutputName '.png']);
imwrite(OutputImage,[OutputName '.jpg'],'Quality',95);

end
