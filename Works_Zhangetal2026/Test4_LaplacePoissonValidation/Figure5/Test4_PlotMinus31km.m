%%
% Test4 plot: TLI line-integration results at -31 km below the surface
% (outside, below the prism).
% Four rows (a,b,c,d), no Fukushima data.  Pan reference is used for the
% absolute-error column.
%%

clc;
clear;
close all;

%%

RE=6371.2;
r=RE-31;
R2=RE;
R1=RE-30;

Fai0=40;
LengthFai=20;
Lamda0=120;
LengthLamda=20;

Fai1=Fai0-LengthFai/2;
Fai2=Fai0+LengthFai/2;
Lamda1=Lamda0-LengthLamda/2;
Lamda2=Lamda0+LengthLamda/2;

%% read TLI results (20 columns, Work1 layout)

fid=fopen('../ResultsPan_-31km.dat','rb+');

NumCoordinates=fread(fid,1,'int32');

Longitude_Pan=fread(fid,NumCoordinates,'double');
Latitude_Pan=fread(fid,NumCoordinates,'double');
Radius_Pan=fread(fid,NumCoordinates,'double');

V_Pan=fread(fid,NumCoordinates,'double');
Vx_Pan=fread(fid,NumCoordinates,'double');
Vy_Pan=fread(fid,NumCoordinates,'double');
Vz_Pan=fread(fid,NumCoordinates,'double');
Vxx_Pan=fread(fid,NumCoordinates,'double');
Vxy_Pan=fread(fid,NumCoordinates,'double');
Vyy_Pan=fread(fid,NumCoordinates,'double');
Vzx_Pan=fread(fid,NumCoordinates,'double');
Vzy_Pan=fread(fid,NumCoordinates,'double');
Vzz_Pan=fread(fid,NumCoordinates,'double');
Vxxx_Pan=fread(fid,NumCoordinates,'double');
Vxxy_Pan=fread(fid,NumCoordinates,'double');
Vxxz_Pan=fread(fid,NumCoordinates,'double');
Vxyz_Pan=fread(fid,NumCoordinates,'double');
Vyyx_Pan=fread(fid,NumCoordinates,'double');
Vyyy_Pan=fread(fid,NumCoordinates,'double');
Vyyz_Pan=fread(fid,NumCoordinates,'double');
Vzzx_Pan=fread(fid,NumCoordinates,'double');
Vzzy_Pan=fread(fid,NumCoordinates,'double');
Vzzz_Pan=fread(fid,NumCoordinates,'double');

fclose(fid);


%% reshape

Longitude_Pan=reshape(Longitude_Pan,401,401);
Latitude_Pan=reshape(Latitude_Pan,401,401);

Vxx_Pan=reshape(Vxx_Pan,401,401);
Vyy_Pan=reshape(Vyy_Pan,401,401);
Vzz_Pan=reshape(Vzz_Pan,401,401);
Vxxx_Pan=reshape(Vxxx_Pan,401,401);
Vyyx_Pan=reshape(Vyyx_Pan,401,401);
Vzzx_Pan=reshape(Vzzx_Pan,401,401);
Vxxy_Pan=reshape(Vxxy_Pan,401,401);
Vyyy_Pan=reshape(Vyyy_Pan,401,401);
Vzzy_Pan=reshape(Vzzy_Pan,401,401);
Vxxz_Pan=reshape(Vxxz_Pan,401,401);
Vyyz_Pan=reshape(Vyyz_Pan,401,401);
Vzzz_Pan=reshape(Vzzz_Pan,401,401);

%% trace quantities (theoretical: 0 outside source, -4*pi inside)

Laplace2nd=Vxx_Pan+Vyy_Pan+Vzz_Pan;
Laplace3rdX=Vxxx_Pan+Vyyx_Pan+Vzzx_Pan;
Laplace3rdY=Vxxy_Pan+Vyyy_Pan+Vzzy_Pan;
Laplace3rd=Vxxz_Pan+Vyyz_Pan+Vzzz_Pan;

Laplace2ndTheory=zeros(401,401);
if ((r>R1)&&(r<R2))
    Index=(Longitude_Pan>Lamda1).*(Longitude_Pan<Lamda2).*(Latitude_Pan>Fai1).*(Latitude_Pan<Fai2);
    Index=logical(Index);
    Laplace2ndTheory(Index)=-4*pi;
end

Laplace2ndError=log10(abs(Laplace2ndTheory-Laplace2nd));
Laplace3rdXError=log10(abs(Laplace3rdX));
Laplace3rdYError=log10(abs(Laplace3rdY));
Laplace3rdError=log10(abs(Laplace3rd));

%% plot data (4 rows, 5 columns)

PlotData={ ...
    Vxx_Pan,Vyy_Pan,Vzz_Pan,Laplace2nd,Laplace2ndError; ...
    Vxxx_Pan,Vyyx_Pan,Vzzx_Pan,Laplace3rdX,Laplace3rdXError; ...
    Vxxy_Pan,Vyyy_Pan,Vzzy_Pan,Laplace3rdY,Laplace3rdYError; ...
    Vxxz_Pan,Vyyz_Pan,Vzzz_Pan,Laplace3rd,Laplace3rdError};

PlotTest4FourRows(Longitude_Pan,Latitude_Pan,PlotData, ...
    [-3,+3],[-0.05,+0.05],[-0.1,+0.1],[-0.5,+0.5], ...
    'Test4_-31km');

return;
