
clc; clear;
mu = 0.5;                         % Mean photon number per pulse
alpha = 0.2;                      % Fiber attenuation (dB/km)
L = 0:0.5:250;                    % Fiber length (km)
t_B = 0.5;                        % Beam splitter transmittance
loss_MZI = 1 * 10^(-0.2);         % 2 dB MZI loss
theta = 0;                        % Phase difference
eta = 0.2;                        % Detector efficiency
P_d = 2e-7;                       % Dark count probability per pulse
t_d = 10e-9;                      % Detector dead time (s)
t_BSE = 0.8;                      % Eve's beam splitter transmittance
f_decoy = 0.2;                    % Decoy-state probability
duty_cycle_TTL = 0.75;            % TTL duty cycle
f = 1.16;                         % Error correction inefficiency
delta_f_values = [1e3 1e6];
R_b_values = [1e8 1e9];
color_ind = {[1 0 0], [1 0.5 0], [0.6 0 0.8],  [0.5 0.25 0] };
color_col = { [0 0 1], [0 0.6 0],[0 0.75 0.75],[0.8 0 0.5] };
labels = {
    '\Delta f=10^3 Hz, R_b=10^8'
    '\Delta f=10^6 Hz, R_b=10^8'
    '\Delta f=10^3 Hz, R_b=10^9'
    '\Delta f=10^6 Hz, R_b=10^9'
};
H2 = @(x) -(x .* log2(max(x,eps)) + (1-x) .* log2(max(1-x,eps)));
figure;
hold on;
h_ind = gobjects(1,4);
h_col = gobjects(1,4);
curve_number = 0;
for k = 1:length(R_b_values)
    R_b = R_b_values(k);
    for j = 1:length(delta_f_values)
    delta_f = delta_f_values(j);
    curve_number = curve_number + 1;
% Phase Noise at LASER
sigma = sqrt(2*pi*delta_f.*(1./R_b));
%Channel transmittance
eta_channel = 10.^(- alpha * L /10);
Photon_count_after_channel = mu * 10.^(- alpha * L /10) .* R_b;
%% DATA LINE %%
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10).* R_b * t_B;
% Photon Counts considering the SPD efficiency and dark count rate
Photon_count_after_D_d = Photon_count_after_BS .* eta + P_d .* R_b;
% Effective time difference between two consecutive pulses at data line
effective_time_difference_D_d = 1./(Photon_count_after_D_d) + t_d;
% Considering the Duty Cycle of TTL Pulses at D_d
effective_counts_D_d = ((1./effective_time_difference_D_d)).* duty_cycle_TTL;
% QBER
Q_data_0 = (P_d*R_b - P_d*R_b*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d);
Q_data = (P_d*R_b - P_d*R_b*(1- exp(-mu.*eta_channel*eta*t_B)))./(effective_counts_D_d).*100;
%% MONITORING LINE%%
% After Beam Splitter
Photon_count_after_BS = mu * 10.^(- alpha * L /10).* R_b * (1-t_B);
% Photon Counts after MZI
Photon_Count_MZI_Upper = loss_MZI * mu * 10.^(- alpha * L /10).*  R_b * (1-t_B) * ((1+ 1*exp(-(sigma.^2)./2) * cos(theta))/2);
Photon_Count_MZI_Lower = loss_MZI * mu * 10.^(- alpha * L /10).*  R_b * (1-t_B) * ((1- 1*exp(-(sigma.^2)./2) * cos(theta))/2);
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
% No of counts after Sifting
R_s = effective_counts_D_d * (1 - f_decoy);
% Eve's Information
I_Eve_Individual =  mu .* 10.^(- 0.2*1/10) .* (1 -  t_BSE) + (1 - V) .*(1 + exp(-mu .* eta_channel)) ./(2*exp(-mu .* eta_channel));
%Overlap
Delta = (2.*V - 1).*exp(-mu) - 2.*sqrt(1 - exp(-2.*mu)).*sqrt(V.*(1 - V));
%Holevo Information under restricted cillective attacks
zeta_COW = Q_data_0 + (1 - Q_data_0).*H2((1 + Delta)./2);
% SKR under selected Individual attack
R_ind = R_s .*(1 - f.*H2(Q_data_0) - I_Eve_Individual);
 R_ind(R_ind < 1e-10) = NaN;
% SKR under restricted Collective attacks
 R_col = R_s .*(1 - f.*H2(Q_data_0) - zeta_COW);
 R_col(R_col < 1e-10) = NaN;
 h_ind(curve_number) = semilogy(L, R_ind,'Color', color_ind{curve_number},  'LineStyle', '-', 'LineWidth', 2,  'DisplayName', labels{curve_number});
 h_col(curve_number) = semilogy(L, R_col,  'Color', color_col{curve_number}, 'LineStyle', '-.', 'LineWidth', 2,'DisplayName', labels{curve_number});

    end
end
xlabel('Length of Fiber (Km)',  'FontName','Times New Roman',  'FontSize',12,'FontWeight','bold');
ylabel('Secret Key Rate (bits/s)', 'FontName','Times New Roman','FontWeight','bold');
grid on;
box on;
xlim([10 250]);
ylim([5e-2 1e8]);
ax_main = gca;
box on;
set(gca, 'YScale','log', 'FontName','Times New Roman', 'FontSize',12, 'LineWidth',2, 'FontWeight','bold');
% Invisible handles for headings
h_ind_title = plot(NaN, NaN, 'LineStyle','none','Marker','none', 'HandleVisibility','on');
h_col_title = plot(NaN, NaN,'LineStyle','none', 'Marker','none', 'HandleVisibility','on');
legend_handles = [
    h_ind_title
    h_ind(:)
    h_col_title
    h_col(:)
];
legend_labels = [
    {'Under Selected Individual attacks'}
    labels
    {'Under Restricted Collective attacks'}
    labels];
legend(legend_handles, legend_labels,'Location','southwest',  'Interpreter','tex','FontSize',9);
hold off;
%% Inset
ax_inset = axes('Position',[0.62 0.60 0.30 0.25]);
copyobj(findobj(ax_main,'Type','line'),ax_inset);
set(ax_inset,'YScale','log');
xlim(ax_inset,[219 222]);
ylim(ax_inset,[1e-1 1e2]);
grid(ax_inset,'on');
box(ax_inset,'on');
set(gca,'FontName','Times New Roman', 'FontSize', 12,'LineWidth', 2, 'FontWeight', 'bold');
hold off;