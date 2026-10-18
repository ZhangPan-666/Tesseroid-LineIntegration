%% Locate Laplace/Poisson trace-error hotspots in TLI results.
clc; clear;
Heights={'1km','-31km','-15km'};
Fai1=30; Fai2=50; La1=110; La2=130;
for k=1:3
 h=Heights{k};
 fid=fopen(['ResultsPan_' h '.dat'],'rb');
 Num=fread(fid,1,'int32');
 lon=fread(fid,Num,'double'); lat=fread(fid,Num,'double'); r=fread(fid,Num,'double');
 for c=1:20, F{c}=fread(fid,Num,'double'); end
 fclose(fid);
 L3x=F{11}+F{15}+F{18}; L3y=F{12}+F{16}+F{19}; L3z=F{13}+F{17}+F{20};
 L3=max([abs(L3x) abs(L3y) abs(L3z)],[],2);
 inside=lon>La1 & lon<La2 & lat>Fai1 & lat<Fai2;
 de=min([min(abs(lon-La1),abs(lon-La2)) min(abs(lat-Fai1),abs(lat-Fai2))],[],2);
 fprintf('==== height %s ====\n',h);
 for t=[1e-8 1e-10 1e-12 1e-13]
  bad=L3>t;
  if ~any(bad), fprintf(' |L3|>%.0e: none\n',t); continue; end
  fprintf(' |L3|>%.0e: n=%d medEdge=%.3f inside=%.2f near0.5=%.2f max=%.2e\n', ...
   t,sum(bad),median(de(bad)),mean(inside(bad)),mean(de(bad)<=0.5),max(L3));
 end
 % second-order trace: theory is -4*pi ONLY for the 15km case whose
 % observation radius lies inside the tesseroid radial interval.
 L2=F{5}+F{7}+F{10};
 L2th=zeros(Num,1);
 if strcmp(h,'-15km'), L2th(inside)=-4*pi; end
 L2e=abs(L2-L2th);
 fprintf(' |L2-theory|>1e-8: n=%d max=%.2e\n',sum(L2e>1e-8),max(L2e));
 for t=[1e-10 1e-11 1e-12]
  if ~any(L2e>t), fprintf('  |L2|>%.0e: none\n',t); continue; end
  fprintf('  |L2|>%.0e: n=%d medEdge=%.3f inside=%.2f max=%.2e\n', ...
   t,sum(L2e>t),median(de(L2e>t)),mean(inside(L2e>t)),max(L2e));
 end
 L3in=max(L3(inside)); L3out=max(L3(~inside));
 fprintf(' max|L3| inside=%.2e outside=%.2e\n',L3in,L3out);
end
fprintf('DONE\n');