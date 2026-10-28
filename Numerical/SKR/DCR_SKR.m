clc; clear;
% parameters
mu = 0.5; % Mean Number of Photons/Pulse
alpha = 0.2; % Attenuation in fiber in dB/Km
delta_f = 1e3; % Linewidth in Hz
d_t = 1e-9; % 1/Bit rate
R_b = 1e9; % Bit Rate
L = 0:0.5:300; % Length of Fiber in Km
t_B = 0.5; % Transmittance of BS
loss_MZI = 1 * 10^(-0.2); % 2 dB loss in MZI
theta = 0; % Phase difference of two consecutive pulses
eta = 0.2; % Efficiency of Detector
t_d = 10e-9; % Dead Time in Detector
t_BSE = 0.8; % Transmittance of Eve's BS
f_decoy = 0.2; % Probability used to generate the decoy state
duty_cycle_TTL = 0.75; % Duty cycle of TTL pulses
f = 1.16;
% Different Dark Count Probability/Pulse 
P_d_values = [1e-6 1e-7 1e-8 ];
lines = {'-','--','-.'};
H2 = @(x) -(x .* log2(max(x,eps)) + (1-x) .* log2(max(1-x,eps)));
nPd = length(P_d_values);
nL = length(L);
R_ind_all = NaN(nPd,nL);
R_col_all = NaN(nPd,nL);
for k = 1:nPd
    P_d = P_d_values(k);
% Phase Noise at LASER
sigma = sqrt(2*pi*delta_f*d_t);
%Channel transmittance
eta_channel = 10.^(- alpha * L /10);
Photon_count_after_channel = mu * 10.^(- alpha * L /10) * R_b;
%% DATA LINE %%
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10)* R_b * t_B;
% Photon Counts considering the SPD efficiency and dark count rate
Photon_count_after_D_d = Photon_count_after_BS .* eta + P_d .* R_b;
% Effective time difference between two consecutive pulses at data line
effective_time_difference_D_d = 1./(Photon_count_after_D_d) + t_d;
% Considering the Duty Cycle of TTL Pulses at D_d
effective_counts_D_d = ((1./effective_time_difference_D_d)).* duty_cycle_TTL;
% QBER
Q_data_0 = (P_d*R_b - P_d*R_b*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d);
Q_data = (P_d*R_b - P_d*R_b*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d).*100;
Q_data_0 = min(max(Q_data_0,0),1);
%% MONITORING LINE%%
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10)* R_b * (1-t_B);
% Photon Counts after MZI
Photon_Count_MZI_Upper = loss_MZI * mu * 10.^(- alpha * L /10)*  R_b * (1-t_B) * ((1+ 1*exp(-(sigma.^2)./2) * cos(theta))/2);
Photon_Count_MZI_Lower = loss_MZI * mu * 10.^(- alpha * L /10)*  R_b * (1-t_B) * ((1- 1*exp(-(sigma.^2)./2) * cos(theta))/2);
% No. of Photons reach at the SPDs
Photon_Count_before_SPD_0 = Photon_Count_MZI_Upper;
Photon_Count_before_SPD_1 = Photon_Count_MZI_Lower;
% Photon Couts considering the SPD efficiency and dark count rate
Photon_Count_after_SPD_0 = Photon_Count_before_SPD_0 * eta + P_d .* R_b;
Photon_Count_after_SPD_1 = Photon_Count_before_SPD_1 * eta + P_d .* R_b;
% Effective time difference between two consecutive pulses
effective_time_difference_SPD_0 = 1./(Photon_Count_after_SPD_0) + t_d;
effective_time_difference_SPD_1 = 1./(Photon_Count_after_SPD_1) + t_d;
% Considering the Duty Cycle of TTL Pulses
effective_counts_SPD_0 = ((1./effective_time_difference_SPD_0))* duty_cycle_TTL;
effective_counts_SPD_1 = ((1./effective_time_difference_SPD_1))* duty_cycle_TTL;
% Considering the Probability factor
total_counts_D_M1 = ((1./effective_time_difference_SPD_0))* duty_cycle_TTL * (2.*f_decoy + (1-f_decoy).^2./2);
total_counts_D_M2 = ((1./effective_time_difference_SPD_1))* duty_cycle_TTL * (2.*f_decoy + (1-f_decoy).^2./2);
% Visibility
V = (total_counts_D_M1 - total_counts_D_M2)./(total_counts_D_M1 + total_counts_D_M2);
V = min(max(V,0),1);
% No of counts after Sifting
R_s = effective_counts_D_d * (1 - f_decoy);
% Eve's Information
I_Eve =  mu .* 10.^(- 0.2*1/10) .* (1 -  t_BSE) + (1 - V) .*(1 + exp(-mu .* eta_channel)) ./(2*exp(-mu .* eta_channel));
%Overlap
Delta = (2.*V - 1).*exp(-mu) - 2.*sqrt(1 - exp(-2.*mu)).*sqrt(V.*(1 - V));
Delta = min(max(Delta,-1),1);
%Holevo Information under restricted cillective attacks
zeta_COW = Q_data_0 + (1 - Q_data_0).*H2((1 + Delta)./2);
% Secret key rate under selected Individual attacks
R_ind = R_s .* (1 - f .* H2(Q_data_0) - I_Eve );
R_ind(R_ind<=0) = eps;
    R_ind_all(k,:) = R_ind;
% Secret key rate under restricted Collective attacks
R_col = R_s .* (1 - f .* H2(Q_data_0) - zeta_COW);
R_col(R_col<=0) = eps;
    R_col_all(k,:) = R_col; 
end
figure(1);
clf;
hold on;
h1 = semilogy(L,R_ind_all(1,:),'LineStyle','-','Color',[0 0.4470 0.7410],  'LineWidth',2);
h2 = semilogy(L,R_ind_all(2,:), 'LineStyle','--','Color',[0 0.4470 0.7410],'LineWidth',2);
h3 = semilogy(L,R_ind_all(3,:),'LineStyle','-.','Color',[0 0.4470 0.7410],  'LineWidth',2);

h4 = semilogy(L,R_col_all(1,:),'LineStyle','-','Color',[0.8500 0.3250 0.0980], 'LineWidth',2);
h5 = semilogy(L,R_col_all(2,:), 'LineStyle','--','Color',[0.8500 0.3250 0.0980], 'LineWidth',2);
h6 = semilogy(L,R_col_all(3,:),'LineStyle','-.','Color',[0.8500 0.3250 0.0980], 'LineWidth',2);
xlabel('Length of Fiber (km)', 'FontName','Times New Roman', 'FontSize',16,'FontWeight','bold');
ylabel('Secret Key Rate (bits/s)', 'FontName','Times New Roman', 'FontSize',16,'FontWeight','bold');
legend([h1 h2 h3 h4 h5 h6], {'P_d = 10^{-6} Hz (Under Selected Individual attack)', ' P_d = 10^{-7} Hz (Under Selected Individual attack)', 'P_d = 10^{-8} Hz (Under Selected Individual attack)', 'P_d = 10^{-6} Hz (Under Restricted Collective attack)', 'P_d = 10^{-7} Hz (Under Restricted Collective attack)',  'P_d = 10^{-8} Hz (Under Restricted Collective attack)'}, 'Location','northeast','FontSize',9);
grid on;
box on;
ax = gca;
set(ax,'YScale','log');
set(ax,'XScale','linear');
xlim([10 300]);
ylim([1e-1 1e8]);
set(gca,'FontName','Times New Roman', 'FontSize',12, 'LineWidth',1.5, 'FontWeight','bold');
hold off;