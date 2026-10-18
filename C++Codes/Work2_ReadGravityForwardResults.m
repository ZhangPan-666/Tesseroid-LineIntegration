% Read tesseroid gravity forward results.

clc;
clear;
close all;

fid=fopen('TFG.ForPar','rb');

AbsTol=fread(fid,1,'double');
RelTol=fread(fid,1,'double');

NumDensityModel=fread(fid,1,'int32');
Fai1=fread(fid,NumDensityModel,'double');
Fai2=fread(fid,NumDensityModel,'double');
Lamda1=fread(fid,NumDensityModel,'double');
Lamda2=fread(fid,NumDensityModel,'double');
R1=fread(fid,NumDensityModel,'double');
R2=fread(fid,NumDensityModel,'double');
rho=fread(fid,NumDensityModel,'double');

fclose(fid);

fid=fopen('Result_TFG.dat','rb');

NumCoordinates=fread(fid,1,'int32');

Longitude=fread(fid,NumCoordinates,'double');
Latitude=fread(fid,NumCoordinates,'double');
Radius=fread(fid,NumCoordinates,'double');

V=fread(fid,NumCoordinates,'double');

Vx=fread(fid,NumCoordinates,'double');
Vy=fread(fid,NumCoordinates,'double');
Vz=fread(fid,NumCoordinates,'double');

Vxx=fread(fid,NumCoordinates,'double');
Vxy=fread(fid,NumCoordinates,'double');
Vyy=fread(fid,NumCoordinates,'double');
Vzx=fread(fid,NumCoordinates,'double');
Vzy=fread(fid,NumCoordinates,'double');
Vzz=fread(fid,NumCoordinates,'double');

Vxxx=fread(fid,NumCoordinates,'double');
Vxxy=fread(fid,NumCoordinates,'double');
Vxxz=fread(fid,NumCoordinates,'double');
Vxyz=fread(fid,NumCoordinates,'double');
Vyyx=fread(fid,NumCoordinates,'double');
Vyyy=fread(fid,NumCoordinates,'double');
Vyyz=fread(fid,NumCoordinates,'double');
Vzzx=fread(fid,NumCoordinates,'double');
Vzzy=fread(fid,NumCoordinates,'double');
Vzzz=fread(fid,NumCoordinates,'double');

fclose(fid);

DataSave=[Longitude,Latitude,Radius,V,Vx,Vy,Vz,...
    Vxx,Vxy,Vyy,Vzx,Vzy,Vzz,...
    Vxxx,Vxxy,Vxxz,Vxyz,Vyyx,...
    Vyyy,Vyyz,Vzzx,Vzzy,Vzzz];
writematrix(DataSave,'GravityForwardResults.dat');
