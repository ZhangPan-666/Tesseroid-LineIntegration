% Create TFG.ForPar for tesseroid gravity forward modelling.

clc;
clear;
close all;

G=1;

Fai0=39;
LengthFai=10;
Lamda0=104;
LengthLamda=10;

Fai1=Fai0-LengthFai/2;
Fai2=Fai0+LengthFai/2;
Lamda1=Lamda0-LengthLamda/2;
Lamda2=Lamda0+LengthLamda/2;

RE=6371.2;

R1=RE-30;
R2=RE;

rho=1.0;

NumDensityModel=size(Lamda0,1);

AbsTol=1e-16;
RelTol=1e-16;

Longitude=104;
Latitude=39;

NumCoordinates=size(Longitude,1);
Radius=(RE+400).*ones(NumCoordinates,1);

FileName='TFG.ForPar';
fid=fopen(FileName,'wb');

fwrite(fid,AbsTol,'double');
fwrite(fid,RelTol,'double');

fwrite(fid,NumDensityModel,'int32');
fwrite(fid,Fai1,'double');
fwrite(fid,Fai2,'double');
fwrite(fid,Lamda1,'double');
fwrite(fid,Lamda2,'double');
fwrite(fid,R1,'double');
fwrite(fid,R2,'double');
fwrite(fid,rho,'double');

fwrite(fid,NumCoordinates,'int32');
fwrite(fid,Longitude,'double');
fwrite(fid,Latitude,'double');
fwrite(fid,Radius,'double');

fclose(fid);
