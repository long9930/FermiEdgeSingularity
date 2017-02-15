clear all;
 
fid = fopen('eigVh.txt');
% fid = fopen('Vhhgrid.txt');
ppp = textscan(fid,'%f');
fclose(fid);
Ng = int32(length(ppp{1}));
Ve_array = zeros(1,Ng);
Ve_array(:) = ppp{1}(1:Ng)';
plot(Ve_array);