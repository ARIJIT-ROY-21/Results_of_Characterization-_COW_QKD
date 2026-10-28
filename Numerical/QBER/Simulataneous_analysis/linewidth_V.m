
clc;
clear;
close all;
mu = 0.5;                  % Mean photon number/pulse
alpha = 0.2;               % Fiber attenuation (dB/km)
L = 80;                    % Fiber length (km)
t_B = 0.5;                 % Beam splitter transmission
loss_MZI = 10^(-0.2);      % MZI insertion loss (2 dB)
theta = 0;                 % Phase difference (radians)
eta = 0.2;                 % Detector efficiency
P_d = 2e-7;                % Dark count probability/pulse
t_d = 10e-9;               % Detector dead time (s)
f_decoy = 0.2;             % Decoy probability
duty_cycle_TTL = 0.75;     % TTL duty cycle
% Laser linewidth in Hz
delta_f_values = logspace(1,8,200);
% Repetition rate in Hz
R_b_values = linspace(1e9,2.5e9,200);
[Delta_F,R_B] = meshgrid(delta_f_values,R_b_values);
d_t = 1./R_B;
sigma_sq = 2*pi.*Delta_F.*d_t;
phase_factor = exp(-sigma_sq./2).*cos(theta);
eta_channel = 10^(-alpha*L/10);
% Photon rate after beam splitter
Photon_count_after_BS = mu.*eta_channel.*R_B.*(1-t_B);
% Upper MZI output
Photon_Count_MZI_Upper = loss_MZI.*Photon_count_after_BS.* (1+phase_factor)./2;
% Lower MZI output
Photon_Count_MZI_Lower = loss_MZI.*Photon_count_after_BS.* (1-phase_factor)./2;
% Detector efficiency and dark counts
Photon_Count_after_SPD_0 =  Photon_Count_MZI_Upper.*eta + P_d.*R_B;
Photon_Count_after_SPD_1 =  Photon_Count_MZI_Lower.*eta + P_d.*R_B;
% Dead-time correction
effective_counts_SPD_0 = Photon_Count_after_SPD_0 ./ (1+Photon_Count_after_SPD_0.*t_d);
effective_counts_SPD_1 = Photon_Count_after_SPD_1 ./(1+Photon_Count_after_SPD_1.*t_d);
effective_counts_SPD_0 = effective_counts_SPD_0.*duty_cycle_TTL;
effective_counts_SPD_1 = effective_counts_SPD_1.*duty_cycle_TTL;
probability_factor = 2*f_decoy + (1-f_decoy)^2/2;
total_counts_D_M1 =  effective_counts_SPD_0.*probability_factor;
total_counts_D_M2 = effective_counts_SPD_1.*probability_factor;
V = (total_counts_D_M1-total_counts_D_M2)./ (total_counts_D_M1+total_counts_D_M2);
fprintf('Minimum Visibility = %.8f\n',min(V(:)));
fprintf('Maximum Visibility = %.8f\n',max(V(:)));
fig = figure( 'Color','w', 'Units','pixels',  'Position',[100 100 1200 850],'Renderer','opengl');
ax = axes('Parent',fig);
set(ax,'Units','normalized');
set(ax,'Position',[0.14 0.19 0.62 0.65]);
h = surf(ax,Delta_F,R_B./1e9,V);
set(h, 'EdgeColor','none', 'FaceColor','interp');
shading interp;
cmap_base = [
    0.000 0.135 0.305
    0.120 0.220 0.400
    0.245 0.320 0.455
    0.360 0.410 0.480
    0.480 0.500 0.470
    0.610 0.595 0.435
    0.745 0.690 0.355
    0.865 0.785 0.275
    0.995 0.910 0.220
];

cmap = interp1(linspace(0,1,size(cmap_base,1)), cmap_base,linspace(0,1,256));
colormap(ax,cmap);
vmin = min(V(:));
vmax = max(V(:));
if vmax > vmin
    caxis(ax,[vmin vmax]);
else
    caxis(ax,[vmin-0.01 vmax+0.01]);
end
set(ax,'XScale','log');
xlim(ax,[1e1 1e8]);
ylim(ax,[1 2.5]);
zlim(ax,[0.7 1]);
set(ax,'XTick',10.^(1:8));
set(ax,'YTick',1:0.25:2.5);
set(ax,'ZTick',0:0.1:1);
set(ax,'XTickLabel',{'10^1','10^2','10^3','10^4','10^5','10^6','10^7','10^8'});
set(ax,'TickLabelInterpreter','tex');
view(ax,45,25);
pbaspect(ax,[1.5 1 0.85]);
hx = xlabel(ax, 'Laser Linewidth (Hz)');
hy = ylabel(ax,'Repetition Rate (GHz)');
hz = zlabel(ax, 'Visibility (V)');
set([hx hy hz], 'FontName','Times New Roman','FontSize',15, 'FontWeight','bold','Interpreter','tex');
set(hz,'Rotation',90);
set(ax, 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','normal', 'LineWidth',1.1, 'Box','on', 'XColor','k',  'YColor','k','ZColor','k');
grid(ax,'on');
cb = colorbar('peer',ax);
set(cb, 'Units','normalized','Position',[0.83 0.24 0.022 0.53], 'FontName','Times New Roman','FontSize',12, 'LineWidth',1);
ylabel(cb,'Visibility, V', 'FontName','Times New Roman', 'FontSize',14,'FontWeight','bold');
set(ax,'Position',[0.14 0.19 0.62 0.65]);