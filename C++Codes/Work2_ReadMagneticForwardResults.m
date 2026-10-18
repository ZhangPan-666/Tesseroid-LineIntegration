% Read tesseroid magnetic forward results.

clc;
clear;
close all;

RE=6371.2;

fid=fopen('TFM.ForPar','rb');

AbsTol=fread(fid,1,'double');
RelTol=fread(fid,1,'double');

NumDensityModel=fread(fid,1,'int32');
Fai1=fread(fid,NumDensityModel,'double');
Fai2=fread(fid,NumDensityModel,'double');
Lamda1=fread(fid,NumDensityModel,'double');
Lamda2=fread(fid,NumDensityModel,'double');
R1=fread(fid,NumDensityModel,'double');
R2=fread(fid,NumDensityModel,'double');
Mx=fread(fid,NumDensityModel,'double');
My=fread(fid,NumDensityModel,'double');
Mz=fread(fid,NumDensityModel,'double');

NumCoordinates=fread(fid,1,'int32');
Longitude=fread(fid,NumCoordinates,'double');
Latitude=fread(fid,NumCoordinates,'double');
Radius=fread(fid,NumCoordinates,'double');

fclose(fid);

fid=fopen('Result_TFM.dat','rb');

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

fclose(fid);

DataSave=[Longitude,Latitude,Radius,V,Vx,Vy,Vz,...
    Vxx,Vxy,Vyy,Vzx,Vzy,Vzz];
writematrix(DataSave,'MagneticForwardResults.dat');
