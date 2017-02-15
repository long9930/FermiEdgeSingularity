clear all;


structurename = 'FES2D';

n2D = 1e12;
Tl = 5;
Te = 5;

prefix = strcat(structurename,'_n2D',num2str(n2D,'%3.2e'),'_Tl',num2str(Tl,'%.0f'),'_Te',num2str(Te,'%.0f'));

fid = fopen(strcat(prefix,'_q.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
q_array = zeros(1,Nq);
q_array(:) = ppp{1}(1:Nq)';

fid = fopen(strcat(prefix,'_Vsq.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
Vsq_array = zeros(1,Nq);
Vsq_array(:) = ppp{1}(1:Nq)';

fid = fopen(strcat(prefix,'_VsPPq1.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
VsPPq1_array = zeros(1,Nq);
VsPPq1_array(:) = ppp{1}(1:Nq)';

fid = fopen(strcat(prefix,'_VsPPq4.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
VsPPq4_array = zeros(1,Nq);
VsPPq4_array(:) = ppp{1}(1:Nq)';

fid = fopen(strcat(prefix,'_PolEq.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
PolEq_array = zeros(1,Nq);
PolEq_array(:) = ppp{1}(1:Nq)';

fid = fopen(strcat(prefix,'_PolHq.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
PolHq_array = zeros(1,Nq);
PolHq_array(:) = ppp{1}(1:Nq)';


fid = fopen(strcat(prefix,'_VsTFq.txt'));
ppp = textscan(fid,'%f');
fclose(fid);
Nq = int32(length(ppp{1}));
VsTFq_array = zeros(1,Nq);
VsTFq_array(:) = ppp{1}(1:Nq)';

% scrsz = get(0,'ScreenSize');
% h = figure('Position',[1 1 0.3*scrsz(3) 0.5*scrsz(4)],...
%      'PaperPosition', [0.25,2,8,8], ...
%      'PaperSize',[8.5 11]);


% h = figure('Position',[1 1 700 600],...
%      'PaperPosition', [0.25,2,8,8], ...
%      'PaperSize',[8.5 11]);

h = figure('Position',[1 1 700 600]);



Vs0 = Vsq_array(1);
q_array = q_array.*1e-9;

% plot(q_array, Vsq_array./Vs0 ,'b', q_array, VsPPq1_array./Vs0,'r:', q_array, VsPPq4_array./Vs0,'c--', q_array, VsTFq_array./Vs0, 'g-.', 'LineWidth',2.5);
% legend('Lindhard','plasmon-pole (C = 1)', 'plasmon-pole (C = 4)','Thomas-Fermi');

plot(q_array, VsTFq_array./Vs0, 'g-.', q_array, VsPPq1_array./Vs0,'r--', q_array, Vsq_array./Vs0 ,'b',   'LineWidth',2.5);
legend('Thomas-Fermi', 'plasmon-pole','Lindhard');


legend('boxoff');
set(findobj('type','axes'),'fontsize',20);
xlabel('q [10^9 m^{-1}]','FontSize',20);
ylabel('V_s(q) / V_s(0)','fontsize',20);

 print(h, '-depsc2', '-painters', '-r2400', 'screening_comparison.eps') 
% plot(q_array, PolEq_array ,'b--', q_array, PolEq_array + PolHq_array,'r','LineWidth',1.0);
%stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV, k_{max} =', 32, num2str(kmax,'%.1f'), 32,'k_F');
% stitle = strcat('n_{2D} =', 32, num2str(n2D,'%3.2e'), 32, 'cm^{-2},',32, 'T_e =',32, num2str(Te,'%.0f'), 32, 'K,', 32, '\gamma =', 32, num2str(gamma,'%.1f'), 32, 'meV');
%title(stitle);
% grid on;


