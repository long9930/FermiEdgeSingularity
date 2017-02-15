clear all;


% structurename = 'testFESnoMagnetic';
structurename = 'FES';

n2D_array = [1e12];
% n2D_array = [0.1e12];

h = figure('Position',[1 1 700 600]);
% hold on;
for kn = 1:length(n2D_array)
    n2D = n2D_array(kn);
Tl = 5;
Te = 5;
gamma = 2;
% kmax = 15;
kkmax = 3.0*sqrt(2.0*pi*5e12*1e4)/2.05;
kF = sqrt(2.0*pi*n2D*1e4);
kmax = kkmax/kF;

prefix = strcat(structurename,'_n2D',num2str(n2D,'%3.2e'),'_Tl',num2str(Tl,'%.0f'),'_Te',num2str(Te,'%.0f'),'_gamma',num2str(gamma,'%.1f'),'_kmax',num2str(kmax,'%.1f'));

fid = fopen(strcat(prefix,'_GainPH.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nph = int32(length(ppp{1})/2);
Eph = zeros(1,Nph);
Gainph = zeros(1,Nph);
Eph(:) = ppp{1}(1:Nph)';
Gainph(:) = ppp{1}(Nph+1:2*Nph)';

fid = fopen(strcat(prefix,'_GainPH0.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nph = int32(length(ppp{1})/2);
Eph0 = zeros(1,Nph);
Gainph0 = zeros(1,Nph);
Eph0(:) = ppp{1}(1:Nph)';
Gainph0(:) = ppp{1}(Nph+1:2*Nph)';

Lp = (8+15)*1e-7;
overlap = 0.01;
Gainph = Gainph ./ Lp .* overlap;
Gainph0 = Gainph0 ./ Lp .* overlap;
% axis auto;
plot(Eph0,Gainph0,'b--', Eph,Gainph,'r','LineWidth',2.0);
end
% legend('free carrier','with Coulomb');
% legend('boxoff');
set(findobj('type','axes'),'fontsize',15);
xlabel('photon energy (eV)','fontsize',15);
ylabel('gain (cm^{-1})','fontsize',15);
% set(gca,'YTick',-0.02:0.01:0.02)
% stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV, k_{max} =', 32, num2str(kmax,'%.1f'), 32,'k_F');
% stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV');
% title(stitle);
% grid on;

 axis([1.30 1.46 -100 100]);

% hold off;

Fig5data = zeros(Nph,3);
Fig5data(:,1) = Eph0';
Fig5data(:,2) = Gainph0';
Fig5data(:,3) = Gainph';

xlswrite('fig5a.xls', Fig5data);


