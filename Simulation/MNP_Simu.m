clc; clear;
% parameters
mu = 0.1:0.01:1; % Mean Number of Photons/Pulse
alpha = 0.2; % Attenuation in fiber in dB/Km
delta_f = 1e3; % Linewidth in Hz
d_t = 1e-9; % 1/Bit rate
R_b = 262144; % Sequence length
L = 40; % Length of Fiber
t_B = 0.5; % Transmittance of BS
loss_MZI = 1 * 10^(0); % 2 dB loss in MZI
theta = 0; % Phase difference of two consecutive pulses
eta = 0.6; % Efficiency of Detector
P_d_counts = 250; % Dark Counts in Detector
t_d = 1e-9; % Dead Time in Detector
duty_cycle_TTL = 1; % Duty cycle of TTL pulses
% Phase Noise at LASER
sigma = sqrt(2*pi*delta_f*d_t);
%Channel transmittance
eta_channel = 10.^(- alpha * L /10);
Photon_count_after_channel = mu * 10.^(- alpha * L /10) * R_b;
% DATA LINE 
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10)* R_b * t_B;
% Photon Counts considering the SPD efficiency and dark count rate
Photon_count_after_D_d = Photon_count_after_BS .* eta + P_d_counts;
% Effective time difference between two consecutive pulses at data line
effective_time_difference_D_d = 1./(Photon_count_after_D_d) + t_d;
% Considering the Duty Cycle of TTL Pulses at D_d
effective_counts_D_d = ((1./effective_time_difference_D_d)).* duty_cycle_TTL;
% QBER
Q_data_0 = (P_d_counts - P_d_counts*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d);
Q_data = (P_d_counts - P_d_counts*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d).*100;
% MONITORING LINE
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10)* R_b * (1-t_B);
% Photon Counts after MZI
Photon_Count_MZI_Upper = loss_MZI * mu * 10.^(- alpha * L /10)*  R_b * (1-t_B) * ((1+ 1*exp(-(sigma.^2)./2) * cos(theta))/2);
Photon_Count_MZI_Lower = loss_MZI * mu * 10.^(- alpha * L /10)*  R_b * (1-t_B) * ((1- 1*exp(-(sigma.^2)./2) * cos(theta))/2);
% No. of Photons reach at the SPDs
Photon_Count_before_SPD_0 = Photon_Count_MZI_Upper;
Photon_Count_before_SPD_1 = Photon_Count_MZI_Lower;
% Photon Couts considering the SPD efficiency and dark count rate
Photon_Count_after_SPD_0 = Photon_Count_before_SPD_0 * eta + P_d_counts;
Photon_Count_after_SPD_1 = Photon_Count_before_SPD_1 * eta + P_d_counts;
% Effective time difference between two consecutive pulses
effective_time_difference_SPD_0 = 1./(Photon_Count_after_SPD_0) + t_d;
effective_time_difference_SPD_1 = 1./(Photon_Count_after_SPD_1) + t_d;
% Considering the Duty Cycle of TTL Pulses
effective_counts_SPD_0 = ((1./effective_time_difference_SPD_0))* duty_cycle_TTL;
effective_counts_SPD_1 = ((1./effective_time_difference_SPD_1))* duty_cycle_TTL;
% Considering the Probability factor
total_counts_D_M1 = ((1./effective_time_difference_SPD_0))* duty_cycle_TTL;
total_counts_D_M2 = ((1./effective_time_difference_SPD_1))* duty_cycle_TTL;
% Visibility
V = (total_counts_D_M1 - total_counts_D_M2)./(total_counts_D_M1 + total_counts_D_M2);
% Excel data
filename = 'Mean_Number_of_photons_simulated_3.xlsx';
data = readtable(filename);
mu_excel = data{:,2};
Q_excel  = data{:,4};
V_excel = data{:,7};
mu_excel = mu_excel(:);
Q_excel = Q_excel(:);
V_excel = V_excel(:);
figure('Color','w');
yyaxis left
p1 = plot(mu, V, 'b-','LineWidth', 3, 'DisplayName', 'Numerical Model (Visibility)');
hold on;
p2 = plot(mu_excel, V_excel, 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'Simulation Result (Visibility)');
ylabel('Visibility (V)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0.68 1]);
yyaxis right
p3 = plot(mu, Q_data, 'r-', 'LineWidth', 3, 'DisplayName', 'Numerical Model (QBER)');
hold on;
p4 = plot(mu_excel, Q_excel, 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'Simulation Result (QBER)');
ylabel('QBER (%)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0 20]);
xlabel('Mean Number of Photons/Pulse', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
xlim([0.1 1]);
legend([p1 p2 p3 p4], 'Location', 'east');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 2);
grid on;
box on;
hold off;
 
