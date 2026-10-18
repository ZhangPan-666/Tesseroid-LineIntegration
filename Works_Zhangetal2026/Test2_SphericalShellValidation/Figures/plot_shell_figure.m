function plot_shell_figure(dataFile, outputBase)
%PLOT_SHELL_FIGURE  TLI-MMSI spherical-shell error maps.
%
%   plot_shell_figure(dataFile, outputBase)
%
%   Strictly follows the panel style of Figure 4 in Zhang et al. (2025)
%   ("ResultsN_Pan*.png"): one Figure with 20 subplots arranged in a
%   4x5 grid, each showing log10 |scaled error| of one gravitational
%   component against the closed-form spherical-shell solution, with the
%   '(a)'--'(t)' labels and the analytical value of each component.
%
%   dataFile   : manuscript Work-format binary .dat produced by
%                shell_validation.cpp (int32 count + 23 double arrays:
%                lon, lat, r, V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy,
%                Vzz, Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy, Vyyz,
%                Vzzx, Vzzy, Vzzz).
%   outputBase : base name ('Results1_Pan1km') producing
%                <base>.png and <base>.jpg.

% ---------- read the Work-format .dat (lon-major, 360 lon x 180 lat) ----------
fid = fopen(dataFile, 'rb');
if fid < 0
    error('Cannot open %s', dataFile);
end
NumCoordinates = fread(fid, 1, 'int32');
if isempty(NumCoordinates) || NumCoordinates ~= 64800
    fclose(fid);
    error('Expected a complete 360x180 file (64800 points), got %d: %s', ...
        NumCoordinates, dataFile);
end
Longitude = fread(fid, NumCoordinates, 'double');
Latitude  = fread(fid, NumCoordinates, 'double');
Radius    = fread(fid, NumCoordinates, 'double');

V    = fread(fid, NumCoordinates, 'double');
Vx   = fread(fid, NumCoordinates, 'double');
Vy   = fread(fid, NumCoordinates, 'double');
Vz   = fread(fid, NumCoordinates, 'double');
Vxx  = fread(fid, NumCoordinates, 'double');
Vxy  = fread(fid, NumCoordinates, 'double');
Vyy  = fread(fid, NumCoordinates, 'double');
Vzx  = fread(fid, NumCoordinates, 'double');
Vzy  = fread(fid, NumCoordinates, 'double');
Vzz  = fread(fid, NumCoordinates, 'double');
Vxxx = fread(fid, NumCoordinates, 'double');
Vxxy = fread(fid, NumCoordinates, 'double');
Vxxz = fread(fid, NumCoordinates, 'double');
Vxyz = fread(fid, NumCoordinates, 'double');
Vyyx = fread(fid, NumCoordinates, 'double');
Vyyy = fread(fid, NumCoordinates, 'double');
Vyyz = fread(fid, NumCoordinates, 'double');
Vzzx = fread(fid, NumCoordinates, 'double');
Vzzy = fread(fid, NumCoordinates, 'double');
Vzzz = fread(fid, NumCoordinates, 'double');
fclose(fid);

allArrays = {Longitude,Latitude,Radius,V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz,...
    Vxxx,Vxxy,Vxxz,Vxyz,Vyyx,Vyyy,Vyyz,Vzzx,Vzzy,Vzzz};
if any(cellfun(@(x) numel(x) ~= NumCoordinates, allArrays))
    error('Truncated or malformed Work-format file: %s', dataFile);
end
if any(~isfinite(Longitude)) || any(~isfinite(Latitude)) || any(~isfinite(Radius))
    error('Non-finite coordinates in: %s', dataFile);
end

% ---------- reshape to lat(rows) x lon(cols) matrices ----------
% Data are stored lon-major: point index = ilon*nlat + ilat.
% reshape(...,180,360) -> column ilon contains all 180 latitudes,
% row ilat is the latitude line, column ilon is the longitude line.
nlat = 180;
nlon = 360;
LonMatrix = reshape(Longitude, nlat, nlon);   % rows: lat, cols: lon
LatMatrix = reshape(Latitude,  nlat, nlon);
V    = reshape(V,    nlat, nlon);
Vx   = reshape(Vx,   nlat, nlon);
Vy   = reshape(Vy,   nlat, nlon);
Vz   = reshape(Vz,   nlat, nlon);
Vxx  = reshape(Vxx,  nlat, nlon);
Vxy  = reshape(Vxy,  nlat, nlon);
Vyy  = reshape(Vyy,  nlat, nlon);
Vzx  = reshape(Vzx,  nlat, nlon);
Vzy  = reshape(Vzy,  nlat, nlon);
Vzz  = reshape(Vzz,  nlat, nlon);
Vxxx = reshape(Vxxx, nlat, nlon);
Vxxy = reshape(Vxxy, nlat, nlon);
Vxxz = reshape(Vxxz, nlat, nlon);
Vxyz = reshape(Vxyz, nlat, nlon);
Vyyx = reshape(Vyyx, nlat, nlon);
Vyyy = reshape(Vyyy, nlat, nlon);
Vyyz = reshape(Vyyz, nlat, nlon);
Vzzx = reshape(Vzzx, nlat, nlon);
Vzzy = reshape(Vzzy, nlat, nlon);
Vzzz = reshape(Vzzz, nlat, nlon);

% 1-deg grid coordinates (identical for every row/column)
LonAxis = LonMatrix(1, :);          % 360 values, -179.5 .. 179.5
LatAxis = LatMatrix(:, 1);          % 180 values, -89.5 .. 89.5
if max(abs(LonMatrix - LonAxis), [], 'all') > 1e-12 || ...
        max(abs(LatMatrix - LatAxis), [], 'all') > 1e-12 || ...
        any(diff(LonAxis) <= 0) || any(diff(LatAxis) <= 0) || ...
        numel(unique(Radius)) ~= 1
    error('Coordinate layout is not a complete lon-major 360x180 grid: %s', ...
        dataFile);
end

% ---------- spherical-shell analytical reference (Appendix B) ----------
RE = 6371.2;
R2 = RE;
R1 = RE - 30;
r  = Radius(1);

if r >= R2
    M0 = 4*pi/3*(R2-R1)*(R2^2+R2*R1+R1^2);
    VSphericalShell    = M0/r;
    VzSphericalShell   = -M0/r^2;
    VxxSphericalShell  = -M0/r^3;
    VyySphericalShell  = -M0/r^3;
    VzzSphericalShell  = 2*M0/r^3;
    VxxzSphericalShell = 3*M0/r^4;
    VyyzSphericalShell = 3*M0/r^4;
    VzzzSphericalShell = -6*M0/r^4;
elseif r <= R1
    VSphericalShell    = 2*pi*(R2-R1)*(R2+R1);
    VzSphericalShell   = 0;
    VxxSphericalShell  = 0;
    VyySphericalShell  = 0;
    VzzSphericalShell  = 0;
    VxxzSphericalShell = 0;
    VyyzSphericalShell = 0;
    VzzzSphericalShell = -0;
else
    h = r - R1;
    radialDifference = h*(r^2+r*R1+R1^2);
    VSphericalShell    = 2*pi*((R2-R1)*(R2+R1)-h^2*(r+2*R1)/(3*r));
    VzSphericalShell   = -4*pi/3*radialDifference/r^2;
    VxxSphericalShell  = -4*pi/3*radialDifference/r^3;
    VyySphericalShell  = VxxSphericalShell;
    VzzSphericalShell  = -4*pi/3*(1+2*R1^3/r^3);
    VxxzSphericalShell = -4*pi*R1^3/r^4;
    VyyzSphericalShell = -4*pi*R1^3/r^4;
    VzzzSphericalShell = 8*pi*R1^3/r^4;
end

VxSphericalShell  = 0;
VySphericalShell  = 0;
VxySphericalShell = 0;
VzxSphericalShell = 0;
VzySphericalShell = 0;
VxxxSphericalShell = 0;
VxxySphericalShell = 0;
VxyzSphericalShell = 0;
VyyxSphericalShell = 0;
VyyySphericalShell = 0;
VzzxSphericalShell = 0;
VzzySphericalShell = 0;

% ---------- absolute errors (log10) ----------
V_Error    = log10(abs(V    - VSphericalShell));
Vx_Error   = log10(abs(Vx   - VxSphericalShell));
Vy_Error   = log10(abs(Vy   - VySphericalShell));
Vz_Error   = log10(abs(Vz   - VzSphericalShell));
Vxx_Error  = log10(abs(Vxx  - VxxSphericalShell));
Vxy_Error  = log10(abs(Vxy  - VxySphericalShell));
Vyy_Error  = log10(abs(Vyy  - VyySphericalShell));
Vzx_Error  = log10(abs(Vzx  - VzxSphericalShell));
Vzy_Error  = log10(abs(Vzy  - VzySphericalShell));
Vzz_Error  = log10(abs(Vzz  - VzzSphericalShell));
Vxxx_Error = log10(abs(Vxxx - VxxxSphericalShell));
Vxxy_Error = log10(abs(Vxxy - VxxySphericalShell));
Vxxz_Error = log10(abs(Vxxz - VxxzSphericalShell));
Vxyz_Error = log10(abs(Vxyz - VxyzSphericalShell));
Vyyx_Error = log10(abs(Vyyx - VyyxSphericalShell));
Vyyy_Error = log10(abs(Vyyy - VyyySphericalShell));
Vyyz_Error = log10(abs(Vyyz - VyyzSphericalShell));
Vzzx_Error = log10(abs(Vzzx - VzzxSphericalShell));
Vzzy_Error = log10(abs(Vzzy - VzzySphericalShell));
Vzzz_Error = log10(abs(Vzzz - VzzzSphericalShell));

% ---------- color bar (identical to manuscript Figure 4) ----------
cbar = MakeColorBar([-18;-15;-12;-9;-6;-3;0],...
    [45;55;77;152;255;255;255],...
    [45;126;175;78;127;255;255],...
    [45;184;74;163;0;51;255]);

% ---------- figure (fix size even in batch mode) ----------
F = figure('Visible', 'off', 'Position', [0 0 1600 880], ...
           'PaperPositionMode', 'auto', 'Color', [1 1 1]);

% data stack in panel order (a..t)
Panels = {V,Vx,Vy,Vz, Vxx,Vxy,Vyy,Vzx,Vzy,Vzz, ...
          Vxxx,Vxxy,Vxxz,Vxyz,Vyyx,Vyyy,Vyyz,Vzzx,Vzzy,Vzzz};
PanelsError = {V_Error,Vx_Error,Vy_Error,Vz_Error, ...
    Vxx_Error,Vxy_Error,Vyy_Error,Vzx_Error,Vzy_Error,Vzz_Error, ...
    Vxxx_Error,Vxxy_Error,Vxxz_Error,Vxyz_Error,Vyyx_Error, ...
    Vyyy_Error,Vyyz_Error,Vzzx_Error,Vzzy_Error,Vzzz_Error};
PanelNames = {'V','Vx','Vy','Vz','Vxx','Vxy','Vyy','Vzx','Vzy','Vzz', ...
    'Vxxx','Vxxy','Vxxz','Vxyz','Vyyx','Vyyy','Vyyz','Vzzx','Vzzy','Vzzz'};
PanelUnits = { '[m^{2}s^{-2}]', '[ms^{-2}]', '[ms^{-2}]', '[ms^{-2}]', ...
    '[s^{-2}]', '[s^{-2}]', '[s^{-2}]', '[s^{-2}]', '[s^{-2}]', '[s^{-2}]', ...
    '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', ...
    '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', ...
    '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]', '[m^{-1}s^{-2}]'};
PanelValues = [VSphericalShell, VxSphericalShell, VySphericalShell, ...
    VzSphericalShell, VxxSphericalShell, VxySphericalShell, ...
    VyySphericalShell, VzxSphericalShell, VzySphericalShell, ...
    VzzSphericalShell, VxxxSphericalShell, VxxySphericalShell, ...
    VxxzSphericalShell, VxyzSphericalShell, VyyxSphericalShell, ...
    VyyySphericalShell, VyyzSphericalShell, VzzxSphericalShell, ...
    VzzySphericalShell, VzzzSphericalShell];
Letters = 'abcdefghijklmnopqrst';

% panel grid positions: 4 rows x 5 columns (paper layout)
RowsPos = [0.77 0.54 0.31 0.08];
ColsPos = [0.05 0.23 0.41 0.59 0.77];

% imagesc coordinate ranges: X = longitude (cols), Y = latitude (rows)
XRange = [LonAxis(1) LonAxis(end)];
YRange = [LatAxis(1) LatAxis(end)];

for k = 1:20
    row = ceil(k/5);
    col = mod(k-1, 5) + 1;

    AxisP = axes('Parent', F);
    AxisP.Position = [ColsPos(col) RowsPos(row) 0.16 0.18];

    IP = imagesc(AxisP, XRange, YRange, PanelsError{k});
    IP.CData = GetCData(PanelsError{k}, cbar, [-18, 0]);
    AxisP.YDir = 'normal';   % latitude increases upward

    AxisP.FontName = 'times new roman';
    AxisP.FontSize = 15;
    AxisP.FontSizeMode = 'manual';
    AxisP.LineWidth = 1.5;
    AxisP.TickDir = 'out';
    AxisP.TickDirMode = 'manual';
    AxisP.XLim = XRange;
    AxisP.YLim = YRange;

    if row == 4
        AxisP.XLabel.String = 'Longitude (^o)';
    else
        AxisP.XTickLabel = [];
    end
    if col == 1
        AxisP.YLabel.String = 'Latitude (^o)';
    else
        AxisP.YTickLabel = [];
    end

    TP = text(AxisP, 0, 0, ['(' Letters(k) ')']);
    TP.FontName = 'times new roman';
    TP.FontWeight = 'normal';
    TP.FontSize = 20;
    TP.Color = [0 0 0];
    TP.Position = [-1.2558*180  1.2007*90  0.0000];

    if PanelValues(k) == 0
        ValueString = sprintf('0 %s', PanelUnits{k});
    else
        Exponent = floor(log10(abs(PanelValues(k))));
        Mantissa = PanelValues(k)/10^Exponent;
        ValueString = sprintf('%.2f\\times10^{%d} %s', Mantissa, Exponent, PanelUnits{k});
    end
    TP = text(AxisP, 0, 0, ['\it{' PanelNames{k} '}\rm=' ValueString]);
    TP.FontName = 'times new roman';
    TP.FontSize = 16;
    TP.Color = [0 0 0];
    TP.Position = [-0.95*180  1.2507*90  0.0000];
end

AxisP.CLim = [-18, 0];

c = colorbar('peer', AxisP);
colormap(AxisP, cbar);
c.Position = [0.9383 0.0771 0.015 0.875];
c.Limits = [-18, 0];
c.Ticks = -18:3:0;
c.TicksMode = 'manual';
c.Label.String = 'log_{10} |scaled error|';
c.Label.FontName = 'times new roman';
c.Label.FontSize = 20;

% ---------- export ----------
print(F, [outputBase '.png'], '-dpng', '-r600');
print(F, [outputBase '.jpg'], '-djpeg95', '-r600');

fprintf('Wrote %s.png / %s.jpg\n', outputBase, outputBase);

end
