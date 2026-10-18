%% Reproduce Figure 2: direct dimensional-reduction verification
% The script reads the two summary_current.csv files in this test and writes
% the PDF/PNG figure into the current figures directory.

close all; clear; clc;
figDir = fileparts(mfilename('fullpath'));
testDir = fileparts(figDir);
cases = { ...
    fullfile(testDir,'results','plus100km','summary_current.csv'), ...
    fullfile(testDir,'results','plus400km','summary_current.csv')};
caseNames = {'100 km','400 km'};

% Prefer Computer Modern when it is installed; otherwise use Times New Roman.
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
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

T = cell(1,2);
for k = 1:2
    T{k} = readtable(cases{k},'TextType','string');
end
components = string(T{1}.component);
x = 1:numel(components);
labels = cellstr(componentLabels(components));
colors = [53 104 168; 192 77 59]/255;
markers = {'o','s'};

fig = figure('Color','w','Units','centimeters','Position',[2 2 18.2 20.0]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
metricNames = {'rel3d2d','rel3d1d','rel3dCpp'};
panelNames = {'(a) 2-D versus 3-D variable precision', ...
              '(b) 1-D versus 3-D variable precision', ...
              '(c) double-precision TLI versus 3-D variable precision'};

for p = 1:3
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    panelValues = [];
    for k = 1:2
        y = double(T{k}.(metricNames{p}));
        y(~isfinite(y) | y <= 0) = realmin('double');
        panelValues = [panelValues; y(:)]; %#ok<AGROW>
        semilogy(ax,x,y,'-','Color',colors(k,:), 'LineWidth',1.0, ...
            'Marker',markers{k},'MarkerSize',4.0,'DisplayName',caseNames{k});
    end
    % Set the logarithmic scale explicitly: do not rely only on semilogy.
    ax.YScale = 'log';
    ax.YLim = decadeLimits(panelValues);
    ax.YLabel.String = 'relative difference';
    ax.FontName = plotFont; ax.FontSize = 9; ax.LineWidth = 0.7;
    ax.TickLabelInterpreter = 'latex';
    ax.XLim = [0.5 numel(x)+0.5];
    ax.YMinorGrid = 'on'; ax.YGrid = 'on'; ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.10;
    text(ax,0.015,0.92,panelNames{p},'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName',plotFont,'FontSize',9,'Interpreter','none', ...
        'BackgroundColor','w','Margin',1);
    legend(ax,'Location','northeast','Box','off','NumColumns',2);
    if p == 3
        ax.XTick = x;
        ax.XTickLabel = labels;
        ax.XTickLabelRotation = 0;
        xlabel(ax,'potential-field quantity');
    else
        ax.XTick = [];
    end
end

% Retain only the two requested raster deliverables.
exportgraphics(fig,fullfile(figDir,'Figure2_DimensionalReduction.png'), ...
    'Resolution',600);
exportgraphics(fig,fullfile(figDir,'Figure2_DimensionalReduction.jpg'), ...
    'Resolution',600);
close(fig);

function labels = componentLabels(components)
labels = strings(size(components));
for i = 1:numel(components)
    c = char(components(i));
    if strcmp(c,'V')
        labels(i) = '$V$';
    else
        labels(i) = ['$V_{' c(2:end) '}$'];
    end
end
end

function lim = decadeLimits(values)
% Return whole-decade limits so differences remain visibly logarithmic.
values = values(isfinite(values) & values > 0);
lo = floor(log10(min(values)));
hi = ceil(log10(max(values)));
if lo == hi
    lo = lo - 1;
    hi = hi + 1;
end
lim = [10^lo 10^hi];
end
