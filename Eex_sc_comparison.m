clear all;

physics_constants;

structurename = 'FES2D';

n2D = 1e12;
Tl = 5;
Te = 5;
kF = sqrt(2*pi*n2D*1e4);

prefix = strcat(structurename,'_n2D',num2str(n2D,'%3.2e'),'_Tl',num2str(Tl,'%.0f'),'_Te',num2str(Te,'%.0f'),'_');

fid = fopen(strcat(prefix,'kk_array.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nkk = int32(length(ppp{1}));
kk_array = zeros(1,Nkk);
kk_array(:) = ppp{1}(1:Nkk)';

fid = fopen(strcat(prefix,'EexLH.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nkk = int32(length(ppp{1}));
EexLH_array = zeros(1,Nkk);
EexLH_array(:) = ppp{1}(1:Nkk)';


fid = fopen(strcat(prefix,'EexPP.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nkk = int32(length(ppp{1}));
EexPP_array = zeros(1,Nkk);
EexPP_array(:) = ppp{1}(1:Nkk)';

fid = fopen(strcat(prefix,'EexTF.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nkk = int32(length(ppp{1}));
EexTF_array = zeros(1,Nkk);
EexTF_array(:) = ppp{1}(1:Nkk)';


h = figure('Position',[1 1 700 600]);

kk_array = kk_array./kF;
plot(kk_array, EexTF_array./echarge*1000, 'g-.', kk_array, EexPP_array./echarge*1000,'r--', kk_array, EexLH_array./echarge*1000, 'b',   'LineWidth',2.5);
legend('Thomas-Fermi', 'plasmon-pole', 'Lindhard');
legend('boxoff');
set(findobj('type','axes'),'fontsize',20);
xlabel('k / k_F','FontSize',20);
ylabel('- E_{ex}(k) [meV]','fontsize',20);

 print(h, '-depsc2', '-painters', '-r2400', 'EexSC_comparison.eps') 
