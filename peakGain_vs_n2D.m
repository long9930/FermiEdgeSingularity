clear all;

structurename = 'FES';

n2D_array = 0.1e12:0.1e12:2.5e12;
Nn2D = length(n2D_array);
EphmaxGain_array = zeros(1,Nn2D);
maxGain_array = zeros(1,Nn2D);


hold on;
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
    
    Lp = (8+15)*1e-7;
    overlap = 0.01;
    Gainph = Gainph ./ Lp .* overlap;
    
    maxGain = 0;
    EphmaxGain = 0;
    for n  = 1:Nph
        if Gainph(n) > maxGain
            maxGain = Gainph(n);
            EphmaxGain = Eph(n);
        end
    end
    
    maxGain_array(kn) = maxGain;
    EphmaxGain_array(kn) = EphmaxGain;


end


h = figure('Position',[1 1 700 600]);

subplot(2,1,1);
plot(n2D_array, maxGain_array,'LineWidth',2.0);
set(findobj('type','axes'),'fontsize',15);
xlabel('electron density (cm^{-2})','fontsize',15);
ylabel('peak gain (cm^{-1})','fontsize',15);

subplot(2,1,2);
plot(n2D_array, EphmaxGain_array,'LineWidth',2.0);
set(findobj('type','axes'),'fontsize',15);
xlabel('electron density (cm^{-2})','fontsize',15);
ylabel('E_{photon} at peak (eV)','fontsize',15);
axis([0 2.5e12 1.33 1.49]);


% ylabel('Gain (cm^{-1})','fontsize',15);
% set(gca,'YTick',-0.02:0.01:0.02)
% stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV, k_{max} =', 32, num2str(kmax,'%.1f'), 32,'k_F');
% stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV');
% title(stitle);
% grid on;


hold off;


Fig5data = zeros(Nn2D,3);
Fig5data(:,1) = n2D_array';
Fig5data(:,2) = maxGain_array';
Fig5data(:,3) = EphmaxGain_array;

xlswrite('fig5b.xls', Fig5data);

