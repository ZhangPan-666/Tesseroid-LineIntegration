% Create TFM.ForPar for tesseroid magnetic forward modelling.

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

NumPermeabilityModel=size(Lamda0,1);

RE=6371.2;

R1=(RE-56).*ones(NumPermeabilityModel,1);
R2=(RE-20).*ones(NumPermeabilityModel,1);
Mx=0.5.*ones(NumPermeabilityModel,1);
My=0.5.*ones(NumPermeabilityModel,1);
Mz=0.5.*ones(NumPermeabilityModel,1);

AbsTol=1e-16;
RelTol=1e-16;

Longitude=104;
Latitude=39;
Radius=RE+1;
NumCoordinates=size(Longitude,1);

FileName='TFM.ForPar';
fid=fopen(FileName,'wb');

fwrite(fid,AbsTol,'double');
fwrite(fid,RelTol,'double');

fwrite(fid,NumPermeabilityModel,'int32');
fwrite(fid,Fai1,'double');
fwrite(fid,Fai2,'double');
fwrite(fid,Lamda1,'double');
fwrite(fid,Lamda2,'double');
fwrite(fid,R1,'double');
fwrite(fid,R2,'double');
fwrite(fid,Mx,'double');
fwrite(fid,My,'double');
fwrite(fid,Mz,'double');

fwrite(fid,NumCoordinates,'int32');
fwrite(fid,Longitude,'double');
fwrite(fid,Latitude,'double');
fwrite(fid,Radius,'double');

fclose(fid);
