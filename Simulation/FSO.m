clc; clear;
% parameters
mu = 0.5074; % Mean Number of Photons/Pulse
alpha_values = [3.8 9.64 28.9]; % Attenuation at Free Space in dB/Km
delta_f = 1e3; % Linewidth in Hz
d_t = 1e-9; % 1/Bit rate
R_b = 262144; % Sequence length
D_t = 0.00005; % Transmitter Aperture in Km
D_r = 0.00005; % Recevier Aperture in Km
D = 0.000025; % Beam Divergence
L_FSO = 0:0.05:2.5; % Transmission Distance
t_B = 0.5; % Transmittance of BS
loss_MZI = 1 * 10^(0); % 2 dB loss in MZI
theta = 0; % Phase difference of two consecutive pulses
eta = 0.6; % Efficiency of Detector
P_d_counts = 250; % Dark Counts in Detector
t_d = 1e-9; % Dead Time in Detector
duty_cycle_TTL = 1; % Duty cycle of TTL pulses
% Phase Noise at LASER
sigma = sqrt(2*pi*delta_f*d_t);
colors = {'r','b',[0.35 0.05 0.45]};
weather_names = {'Haze','Rain','Fog'};
line_styles = {'-','--','-.'};
filename = 'FSO_simu_4.xlsx';
data = readtable(filename);
L_FSO_excel = data{:,1};
Q_excel_haze = data{:,3};
Q_excel_rain  = data{:,5};
Q_excel_fog = data{:,7};
Q_excel_all = {Q_excel_haze, Q_excel_rain, Q_excel_fog};
figure;
hold on;
for k = 1:length(alpha_values)

    alpha = alpha_values(k);
eta_channel = (D_r./(D_t + D.*L_FSO)).^2 .* 10.^(- alpha .* L_FSO /10);
Photon_count_after_channel = mu * eta_channel * R_b;
% DATA LINE 
% After Beam Splitter
Photon_count_after_BS = mu * eta_channel* R_b * t_B;
% Photon Counts considering the SPD efficiency and dark count rate
Photon_count_after_D_d = Photon_count_after_BS .* eta + P_d_counts;
% Effective time difference between two consecutive pulses at data line
effective_time_difference_D_d = 1./(Photon_count_after_D_d) + t_d;
% Considering the Duty Cycle of TTL Pulses at D_d
effective_counts_D_d = ((1./effective_time_difference_D_d)).* duty_cycle_TTL;
% Error Counts at the D_d
E_C_D_d = ((P_d_counts.*exp(-mu.*eta_channel.*eta.*t_B))./(1 + t_d.*P_d_counts.*exp(-mu.*eta_channel.*eta.*t_B)))* duty_cycle_TTL;
% QBER
Q_data_0 = (E_C_D_d) ./(effective_counts_D_d);
Q_data = (E_C_D_d) ./(effective_counts_D_d).*100;
% Plot
plot(L_FSO, Q_data, 'Color', colors{k}, 'LineStyle', line_styles{k}, 'LineWidth', 2.5, 'DisplayName', [weather_names{k} ' - Numerical Model']);
    Q_excel = Q_excel_all{k};
simulation_colors = {'r', 'b', [0.35 0.05 0.45]};
plot(L_FSO_excel, Q_excel,'o','Color', simulation_colors{k},'MarkerFaceColor', simulation_colors{k},'MarkerSize', 7,'LineWidth', 1.5,'DisplayName', [weather_names{k} ' - Simulation Result']);

end
xlabel('Transmission Distance (Km)','FontName', 'Times New Roman', 'FontSize', 13,'FontWeight', 'bold');
ylabel('QBER (%)','FontName', 'Times New Roman', 'FontSize', 13,'FontWeight', 'bold');
legend('Location','northwest');
xlim([0 2.5]);
ylim([0 65]);
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12,'LineWidth', 2);
hold off;