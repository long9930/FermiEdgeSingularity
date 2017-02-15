clear all;

fid1 = fopen('FES_Tl5_Vqfrom.txt');
ppp1 = textscan(fid1,'%f');
fclose(fid1);
Nq1 = int32(length(ppp1{1})/2);
qform_array1 = zeros(1,Nq1);
Vqform_array1 = zeros(1,Nq1);
qform_array1(:) = ppp1{1}(1:Nq1)';
Vqform_array1(:) = ppp1{1}(Nq1+1:2*Nq1)';
plot(qform_array1, Vqform_array1);

% hold on;
% 
% fid = fopen('FES_Tl5_Vqfrom.txt');
% ppp = textscan(fid,'%f');
% fclose(fid);
% Nq = int32(length(ppp{1})/2);
% qform_array = zeros(1,Nq);
% Vqform_array = zeros(1,Nq);
% qform_array(:) = ppp{1}(1:Nq)';
% Vqform_array(:) = ppp{1}(Nq+1:2*Nq)';
% plot(qform_array, Vqform_array);
% 
% plot(qform_array, Vqform_array./Vqform_array1);