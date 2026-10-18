%% Test 3: relative error, absolute error, and value-overlay figures
clc; clear; close all;

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here,'results','polar_validation.mat'));

availableFonts = string(listfonts);
if any(strcmpi(availableFonts,'Computer Modern'))
    plotFont = 'Computer Modern';
elseif any(contains(lower(availableFonts),'computer modern'))
    plotFont = char(availableFonts(find(contains(lower(availableFonts),'computer modern'),1)));
else
    plotFont = 'Times New Roman';
end
set(groot,'defaultAxesFontName',plotFont);
set(groot,'defaultTextFontName',plotFont);
set(groot,'defaultLegendFontName',plotFont);

names = {'V','V_z','V_{xx}','V_{yy}','V_{zz}', ...
         'V_{xxz}','V_{yyz}','V_{zzz}'};
units = {'[m^2 s^{-2}]','[m s^{-2}]','[s^{-2}]','[s^{-2}]', ...
         '[s^{-2}]','[m^{-1} s^{-2}]','[m^{-1} s^{-2}]','[m^{-1} s^{-2}]'};
tags = {'Normal','Small','Thin','ZonalBand'};
black = [0 0 0];
theoryGray = [0.60 0.60 0.60];
computedBlue = [0.00 0.25 0.85];

for k = 1:numel(S.cases)
    makeErrorFigure(S,'relativeError',names,units,tags{k},k, ...
        'RelativeError',[-18 1],[-16 -8 0],black,plotFont);
    makeErrorFigure(S,'absoluteError',names,units,tags{k},k, ...
        'AbsoluteError',[-30 1],[-28 -20 -12 -4 0],black,plotFont);
    makeValueFigure(S,names,units,tags{k},k,theoryGray,computedBlue,plotFont);
end
fprintf('Wrote 12 PNG figures.\n');

function makeErrorFigure(S,field,names,units,tag,k,suffix,xlimValue,xticks,lineColor,plotFont)
    F = figure('Visible','off','Color','white','Position',[1 1 1700 900]);
    tl = tiledlayout(F,1,8,'TileSpacing','compact','Padding','compact');
    for j = 1:8
        ax = nexttile(tl);
        e = log10(max(S.(field)(:,j,k),realmin('double')));
        plot(ax,e,S.height,'Color',lineColor,'LineWidth',1.5);
        formatAxis(ax,names{j},units{j},j,plotFont);
        if strcmp(field,'absoluteError')
            ax.XLim = [-26 -7];
            ax.XTick = [-24 -16 -8];
        else
            ax.XLim = xlimValue;
            ax.XTick = xticks;
        end
        ax.XLabel.String = ['log_{10} |' lower(strrep(suffix,'Error',' error')) '|'];
        ax.XLabel.Interpreter = 'tex';
        if j == 1
            ax.YLabel.String = 'Height relative to $R_E$ [km]';
            ax.YLabel.Interpreter = 'latex';
        end
    end
    sgtitle(tl,sprintf('%s: %s',S.cases(k).name,strrep(suffix,'Error',' error')), ...
        'FontName',plotFont,'FontSize',20,'FontWeight','normal');
    outBase = fullfile(fileparts(mfilename('fullpath')), ...
        ['Test3_' tag '_' suffix]);
        exportgraphics(F,[outBase '.png'],'Resolution',200);
    close(F);
end

function makeValueFigure(S,names,units,tag,k,theoryGray,computedBlue,plotFont)
    F = figure('Visible','off','Color','white','Position',[1 1 1700 900]);
    tl = tiledlayout(F,1,8,'TileSpacing','compact','Padding','compact');
    legendHandles = gobjects(2,1);
    for j = 1:8
        ax = nexttile(tl);
        theory = S.analytical(:,j,k);
        computed = S.numerical(:,j,k);
        hTheory = plot(ax,theory,S.height,'Color',theoryGray,'LineWidth',4);
        hold(ax,'on');
        hComputed = plot(ax,computed,S.height,'Color',computedBlue,'LineWidth',1);
        hold(ax,'off');
        if j == 8
            legendHandles = [hTheory;hComputed];
        end
        formatAxis(ax,names{j},units{j},j,plotFont);
        values = [theory;computed];
        vmin = min(values); vmax = max(values);
        pad = max(1e-30,0.06*max(vmax-vmin,abs(vmax)*1e-12));
        ax.XLim = [vmin-pad,vmax+pad];
        ax.XLabel.String = 'Value';
        ax.XLabel.Interpreter = 'tex';
        if j == 1
            ax.YLabel.String = 'Height relative to $R_E$ [km]';
            ax.YLabel.Interpreter = 'latex';
        end
    end
    % Place the legend in the upper-right corner of the rightmost panel.
    lgd = legend(legendHandles,{'Theory','Computed'},'Location','northeast', ...
        'FontName',plotFont,'FontSize',13,'Box','on');
    sgtitle(tl,sprintf('%s: theory and computed values',S.cases(k).name), ...
        'FontName',plotFont,'FontSize',20,'FontWeight','normal');
    outBase = fullfile(fileparts(mfilename('fullpath')), ...
        ['Test3_' tag '_Values']);
    exportgraphics(F,[outBase '.png'],'Resolution',200);
    close(F);
end

function formatAxis(ax,name,unit,j,plotFont)
    ax.FontName = plotFont;
    ax.FontSize = 14;
    ax.LineWidth = 1.5;
    ax.YLim = [-1000 1000];
    ax.YTick = [-1000 -500 0 500 1000];
    ax.TickDir = 'out';
    ax.Box = 'on';
    ax.XLabel.FontName = plotFont;
    ax.XLabel.FontSize = 16;
    ax.YLabel.FontName = plotFont;
    ax.YLabel.FontSize = 17;
    unitMath = strrep(strrep(unit,'[',''),']','');
    title(ax,['$' name '\;[' unitMath ']$'],'Interpreter','latex', ...
        'FontName',plotFont,'FontSize',18, ...
        'FontWeight','normal','Color','black');
    text(ax,0.04,0.92,sprintf('(%c)',char('a'+j-1)), ...
        'Units','normalized','FontName',plotFont,'FontSize',17, ...
        'FontWeight','normal','Color','black','Interpreter','none');
    if j ~= 1
        ax.YTickLabel = [];
        ax.YLabel.String = '';
    end
end
