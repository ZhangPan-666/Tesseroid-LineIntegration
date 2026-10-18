%% Compare second-order trace residual of TLI vs reference Pan.
% Both methods solve the same tesseroid; if both have the same ~1e-9 trace
% residual near the prism edges, the "poor-looking" value is an intrinsic
% property of the trace cancellation at edges, not a TLI defect.
clc; clear;
Heights={'1km','-31km','-15km'};
Fai1=30; Fai2=50; La1=110; La2=130;
RefDir=getenv('TLI_TRACE_REFERENCE_DIR');
assert(isfolder(RefDir), 'Set TLI_TRACE_REFERENCE_DIR to the external comparison results.');
for k=1:3
 h=Heights{k};
 fid=fopen(['ResultsPan_' h '.dat'],'rb');
 Num=fread(fid,1,'int32');
 lon=fread(fid,Num,'double'); lat=fread(fid,Num,'double');
 r=fread(fid,Num,'double');
 for c=1:20, T{c}=fread(fid,Num,'double'); end
 fclose(fid);
 fid=fopen(fullfile(RefDir,['ResultsPan_' h '.dat']),'rb');
 Num2=fread(fid,1,'int32');
 lon2=fread(fid,Num2,'double'); lat2=fread(fid,Num2,'double');
 r2=fread(fid,Num2,'double');
 for c=1:20, P{c}=fread(fid,Num2,'double'); end
 fclose(fid);
 L2T=T{5}+T{7}+T{10};
 L2P=P{5}+P{7}+P{10};
 % 15km internal: subtract -4pi inside
 if strcmp(h,'-15km')
  in=lon>La1 & lon<La2 & lat>Fai1 & lat<Fai2;
  L2T=L2T+4*pi*in; L2P=L2P+4*pi*in;
 end
 fprintf('==== %s : trace-residual |Vxx+Vyy+Vzz| (theory 0) ====\n',h);
 fprintf('  TLI : max=%.3e  median=%.3e  p99=%.3e\n',max(abs(L2T)),median(abs(L2T)),quantile(abs(L2T),0.99));
 fprintf('  PanRef: max=%.3e  median=%.3e  p99=%.3e\n',max(abs(L2P)),median(abs(L2P)),quantile(abs(L2P),0.99));
end
fprintf('DONE\n');