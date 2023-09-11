%% Iono estimation

%This assumes the variables storeData and model_save have been alreadty
%imported

%% Define various constants

Const;
H = 300e3;
load('/Users/Bernard/Documents/GitHub/raPPPid/RESULTS/COM8_191006_234852_T/noname_28-Apr-2023_09h43m29s/data4plot.mat');

% load color maps
load('/Users/Bernard/Downloads/ScientificColourMaps8/roma/DiscretePalettes/roma10.mat');
load('/Users/Bernard/Downloads/ScientificColourMaps8/bam/DiscretePalettes/bam10.mat');
load('/Users/Bernard/Downloads/ScientificColourMaps8/roma/DiscretePalettes/roma100.mat')


% true if GNSS was processed and should be plotted
isGPS = settings.INPUT.use_GPS;
isGLO = settings.INPUT.use_GLO;
isGAL = settings.INPUT.use_GAL;
isBDS = settings.INPUT.use_BDS;
noGNSS = isGPS + isGLO + isGAL + isBDS;

% satellite indices for each GNSS
idx_gps = 001:000+DEF.SATS_GPS;
idx_glo = 101:100+DEF.SATS_GLO;
idx_gal = 201:200+DEF.SATS_GAL;
% idx_bds = 301:size(bool_obs,2);

gamma_phi_gps = Const.GPS_F2^2 / (Const.GPS_F1^2 - Const.GPS_F2^2);
gamma_rho_gps = Const.GPS_F2^2 / (Const.GPS_F2^2 - Const.GPS_F1^2);

gamma_phi_gal = Const.GAL_F5b^2 / (Const.GAL_F1^2 - Const.GAL_F5b^2);
gamma_rho_gal = Const.GAL_F5b^2 / (Const.GAL_F5b^2 - Const.GAL_F1^2);

gamma_glo = Const.GLO_k2^2 / (Const.GLO_F1^2 - Const.GLO_k2^2);
gamma_bds = Const.BDS_F2^2 / (Const.BDS_F1^2 - Const.BDS_F2^2);

gps_prn = 26;
gal_prn = 33;


elev_thresh = 10;
%% Estimate iono from code and phase
sat_idx = idx_gps(gps_prn);
[cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh);

%% Plot estimates
i_plot = 1;

fig_title = 'Code Iono Estimates: Multiple Sats';
togglefig(fig_title); clf;
gps_sats = [1, 3, 9, 26, 31];
gps_sats = [3, 26, 9]; %[3, 26, 9];
subplot_idx(1) = length(gps_sats);
j=1;
for i = gps_sats
    sat_idx = idx_gps(i);
    [cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh);
    plot_iono_estimates(cmc, cmp, sat_dcb, iono_gim, gamma_rho_gps, gamma_phi_gps, [length(gps_sats), j], sat_idx);
    j = j+1;
end

fig_title = 'Code Iono Estimates vs. Elevation: Multiple Sats';
togglefig(fig_title); clf;
j=1;
for i = gps_sats
    sat_idx = idx_gps(i);
    [cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh);
    plot_iono_vs_elev(cmc, cmp, sat_dcb, iono_gim, gamma_rho_gps, gamma_phi_gps, satellites,...
        [length(gps_sats), j], sat_idx)
    j = j+1;
end

fig_title = 'Receiver DCB';
togglefig(fig_title); clf;
rec_dcb = nanmean(gamma_rho_gps*((1/gamma_rho_gps *iono_gim) - cmc - sat_dcb));
plot(gamma_rho_gps*((1/gamma_rho_gps *iono_gim) - cmc - sat_dcb), 'b.');
hold on;
plot(ones(length(cmc), 1)*rec_dcb, 'r.');
xlabel('Epoch (s)')
ylabel('DCB_r')
grid minor
%%
fig_title = 'Code, Phase Iono Estimates: Multiple Sats';
% fig_c1 = figure('Name', fig_title, 'units','normalized', 'outerposition',[0 0 1 1], 'NumberTitle','off');
togglefig(fig_title); clf;
PLOT = 0;
j = 1;
for i = gps_sats
    sat_idx = idx_gps(i);
    [cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh);
    % hax (j) = subplot(1, length(gps_sats), j);
    hax(j) = subplot(length(gps_sats), 1, j);
    % approx integer ambiguity difference (lambda_1*N_1 - lambda_2*N_2)
    B = nanmean((gamma_phi_gps * cmp) - iono_gim);
    rec_dcb = nanmean(gamma_rho_gps*(cmc + sat_dcb) - iono_gim);
    if PLOT
        hold on;
        grid minor;
        plot(gamma_rho_gps * (cmc + sat_dcb) - rec_dcb, 'Marker', '.', 'Linestyle', '-',  'Color', roma10(9,:));
        plot((gamma_phi_gps * cmp) - B , 'Marker', '.', 'Linestyle', '-',  'Color', [0.9 roma10(1,2:3)]);
        plot(iono_gim, 'Marker', '.', 'Linestyle', '-',  'Color', [bam10(1,1), 0.85, bam10(1,3)]);

        xlabel('Epoch (s)')
        ylabel('Iono Estimate (m)')
        sys = gnss2char('GPS');
        sat_prn = strcat(sys, num2str(mod(sat_idx',100), '%02.0f'));
        title(sat_prn)
        legend('Receiver pseudorange iono estimate', 'Receiver phase iono estimate', 'IGS GIM iono estimate', 'FontSize', 14)
        box on
    end
    j = j+1;
end

linkaxes(hax, 'x')
%%
fig_title = 'Code, Phase VTEC Estimates: Multiple Sats';
% fig_c1 = figure('Name', fig_title, 'units','normalized', 'outerposition',[0 0 1 1], 'NumberTitle','off');
togglefig(fig_title); clf;
j = 1;

plot_x_elev = 1;
PLOT = 0;

vtec_rho = zeros(size(satellites.elev));
vtec_phi = zeros(size(satellites.elev));
vtec_gim = zeros(size(satellites.elev));
gps_sats = 1:32; 
for i = [idx_gps, idx_gal]
    % sat_idx = idx_gps(i);
    sat_idx = i;
    [cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh);
    % approx integer ambiguity difference (lambda_1*N_1 - lambda_2*N_2)

    if isGnssSys(sat_idx, idx_gps)
        gamma_phi = gamma_phi_gps;
        gamma_rho = gamma_rho_gps;
        f1 =  Const.GPS_F1;
    elseif isGnssSys(sat_idx, idx_gal)
        gamma_phi = gamma_phi_gal;
        gamma_rho = gamma_rho_gal;
        f1 =  Const.GAL_F1;
    end

    B = nanmean((gamma_phi * cmp) - iono_gim); % integer ambiguity difference
    rec_dcb = nanmean(gamma_rho*(cmc + sat_dcb) - iono_gim); % receiver DCB

    iono_rho_estimate = gamma_rho * (cmc + sat_dcb) - rec_dcb;
    iono_phi_estimate = (gamma_phi * cmp) - B;

    sat_elev_rad = deg2rad(satellites.elev(:,sat_idx));
    vtec_rho(:,i) = compute_vtec(sat_elev_rad, f1, iono_rho_estimate);
    vtec_phi(:,i) = compute_vtec(sat_elev_rad, f1, iono_phi_estimate);
    vtec_gim(:,i) = compute_vtec(sat_elev_rad, f1, iono_gim);

    if PLOT
        % hax (j) = subplot(length(gps_sats), 1, j);
        hax(j) = subplot(1, length(gps_sats), j);
        hold on;
        grid minor;
        if ~plot_x_elev

            plot(vtec_rho(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', roma10(9,:));
            plot(vtec_phi(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', [0.9 roma10(1,2:3)]);
            plot(vtec_gim(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', [bam10(1,1), 0.85, bam10(1,3)]);
            xlabel('Epoch (s)')
        else
            sat_elev_deg = rad2deg(sat_elev_rad);
            plot(sat_elev_deg, vtec_rho(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', roma10(9,:));
            plot(sat_elev_deg, vtec_phi(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', [0.9 roma10(1,2:3)]);
            plot(sat_elev_deg, vtec_gim(:,i), 'Marker', '.', 'Linestyle', '-',  'Color', [bam10(1,1), 0.85, bam10(1,3)]);
            xlabel('Satellite Elevation (deg)')
            xlim([5, 85]);
        end


        ylabel('VTEC (TECu)')
        ylim([-0, 20])
        sys = gnss2char('GPS');
        sat_prn = strcat(sys, num2str(mod(sat_idx',100), '%02.0f'));
        title(sat_prn)
        legend('Receiver pseudorange VTEC estimate', 'Receiver phase VTEC estimate', 'IGS GIM VTEC', 'FontSize', 14)
        box on
    end
    j = j+1;
end
if PLOT
    linkaxes(hax, 'x');
end

%% Plot the iono pierce points vs elevation 
pos_geo = cart2geo(storeData.param(end,1:3));
[lat_pp, lon_pp] = calcIPP(pos_geo.ph, pos_geo.la, satellites.az*pi/180, satellites.elev*pi/180, H);
%
fig_title = 'Iono PP map: Elevation';
% fig_c1 = figure('Name', fig_title, 'units','normalized', 'outerposition',[0 0 1 1], 'NumberTitle','off');
togglefig(fig_title); clf;
sats_to_plot = 1:32;
NumSamples = 120;
for i = sats_to_plot
    geoscatter(full(lat_pp(1:NumSamples, i))*180/pi, full(lon_pp(1:NumSamples, i))*180/pi, 1, satellites.elev(1:NumSamples,i), '.')
    hold on;
end
geoscatter(pos_geo.ph*180/pi, pos_geo.la*180/pi, 50, 'Marker', 'square', 'MarkerFaceColor','r')
geobasemap darkwater
geolimits([18.0930 50],[-120 -120])
grid off
sat_prn = strcat(sys, num2str(mod(sats_to_plot',100), '%02.0f'));
% legend(sat_prn, 'Fontsize', 16)
h = colorbar;
h.Label.String = "Elevation (deg)";
h.Label.FontSize = 16;

for i = gps_sats
    % geoscatter(full(lat_pp(:, i))*180/pi, full(lon_pp(:, i))*180/pi, 20, satellites.elev(:,i), '+')
    hold on;
end
sat_prn = strcat(sys, num2str(mod(gps_sats',100), '%02.0f'));
colormap(flip(roma100))
geobasemap streets-dark
% legend(sat_prn, 'Fontsize', 16)

rx_elev = (0:20:80)*pi/180;
rx_az = 0:0.01:2*pi;
for elev=rx_elev
    [phi_pp, lam_pp] = calcIPP(pos_geo.ph, pos_geo.la, rx_az, elev, H);
    geoscatter(phi_pp*180/pi, lam_pp*180/pi, 1, ones(1,length(phi_pp))*elev*180/pi, '*')
end

%% lot the iono pierce points vs VTEC 
fig_title = 'Iono PP map: VTEC';
% fig_c1 = figure('Name', fig_title, 'units','normalized', 'outerposition',[0 0 1 1], 'NumberTitle','off');
togglefig(fig_title); clf;
sats_to_plot = [idx_gps, idx_gal];
% geoscatter(pos_geo.ph*180/pi, pos_geo.la*180/pi, 100, 'Marker', 'square', 'MarkerFaceColor','r')
NumSamples = 120;
for i = sats_to_plot
    if isGnssSys(i, idx_gps)
        geoscatter(full(lat_pp(1:NumSamples, i))*180/pi, full(lon_pp(1:NumSamples, i))*180/pi, 100, vtec_phi(1:NumSamples,i), '.');
    elseif isGnssSys(i, idx_gal)
        geoscatter(full(lat_pp(1:NumSamples, i))*180/pi, full(lon_pp(1:NumSamples, i))*180/pi, 100, vtec_phi(1:NumSamples,i), 'x');

    end
    hold on;
end
geoscatter(pos_geo.ph*180/pi, pos_geo.la*180/pi, 100, 'Marker', 'square', 'MarkerFaceColor','r')

geobasemap darkwater
geolimits([18.0930 50],[-120 -120])
grid off
sat_prn = strcat(sys, num2str(mod(sats_to_plot',100), '%02.0f'));
% legend(sat_prn, 'Fontsize', 16)
h = colorbar;
h.Label.String = "VTEC (TECu)";
h.Label.FontSize = 16;

% endpoint = 120;
% for i = gps_sats
%     geoscatter(full(lat_pp(1:endpoint, i))*180/pi, full(lon_pp(1:endpoint, i))*180/pi, 20, vtec_phi(1:endpoint,i), '+')
%     hold on;
% end
sat_prn = strcat(sys, num2str(mod(gps_sats',100), '%02.0f'));
colormap(flip(roma100))
geobasemap streets-dark

% legend(sat_prn, 'Fontsize', 16)
%%
function [cmc, cmp, sat_dcb, iono_gim] = compute_iono_estimates (storeData, satellites, sat_idx, elev_thresh)
%GIM iono
iono_gim = storeData.iono_corr(:, sat_idx);
iono_gim(iono_gim == 0) = NaN;

cmp = (storeData.L1(:, sat_idx) - storeData.L2(:, sat_idx));
cmp = cmp .* (satellites.elev(:, sat_idx)>elev_thresh);
cmp(cmp==0) = NaN;

cmc = (storeData.C1(:, sat_idx) - storeData.C2(:, sat_idx));
% throw out sats below elevation of 10 deg
cmc = cmc .* (satellites.elev(:, sat_idx)>elev_thresh);
cmc(cmc==0) = NaN;

sat_dcb = storeData.C1_bias(:,sat_idx) - storeData.C2_bias(:, sat_idx);
end

function plot_iono_estimates(cmc, cmp, sat_dcb, iono_gim, gamma_rho, gamma_phi, subplot_idx, sat_idx)
subplot(subplot_idx(1), 1, subplot_idx(2))
plot(gamma_rho * (cmc + sat_dcb), 'b.');
hold on;
% plot(gamma_phi * (cmp + sat_dcb), 'r.')
plot(iono_gim, 'r.');
ylim([-15 15])
xlim(flip(size(iono_gim)))
grid minor
xlabel('Epoch (s)')
ylabel('Iono Estimate (m)')
legend('Receiver Pseudorange Iono estimate', 'IGS GIM Iono estimate', 'FontSize', 14)
sys = gnss2char('GPS');
sat_idx = strcat(sys, num2str(mod(sat_idx',100), '%02.0f'));
title(sat_idx)
end

function plot_tec_estimates(cmc, cmp, sat_dcb, iono_gim, gamma_rho, gamma_phi, subplot_idx, sat_idx)
subplot(subplot_idx(1), 1, subplot_idx(2))
plot(gamma_rho * (cmc + sat_dcb), 'b.');
hold on;
% plot(gamma_phi * (cmp + sat_dcb), 'r.')
plot(iono_gim, 'r.');
% ylim([-15 15])
xlim(flip(size(iono_gim)))
grid minor
xlabel('Epoch (s)')
ylabel('Iono Estimate (m)')
legend('Receiver Pseudorange Iono estimate', 'IGS GIM Iono estimate', 'FontSize', 14)
sys = gnss2char('GPS');
sat_idx = strcat(sys, num2str(mod(sat_idx',100), '%02.0f'));
title(sat_idx)
end

function plot_iono_vs_elev(cmc, cmp, sat_dcb, iono_gim, gamma_rho, gamma_phi, satellites, subplot_idx, sat_idx)
subplot(subplot_idx(1), 1, subplot_idx(2));
hold on;
plot(satellites.elev(:, sat_idx), gamma_rho * (cmc + sat_dcb) + 13, 'b.')
plot(satellites.elev(:, sat_idx), iono_gim, 'r.');
grid minor
xlabel('Elevation (deg)', 'FontSize', 14)
ylabel('Iono Estimate (m)', 'FontSize', 14)
legend('Receiver Pseudorange Iono estimate', 'IGS GIM Iono estimate', 'FontSize', 14)
sys = gnss2char('GPS');
sat_idx = strcat(sys, num2str(mod(sat_idx',100), '%02.0f'));
title(sat_idx)
xlim([10 90]);

box on

end

function vtec = compute_vtec(elev, f, iono_estimate)
hI = 350e3; % ionospheric shell height
stec = iono_estimate * f^2/40.3e16;
vtec = stec .* sin(elev + compute_pierce_point_zenith_angle_rad(elev, hI));
% vtec = stec;
end

function psi_pp = compute_pierce_point_zenith_angle_rad(elev, hI)
psi_pp = pi/2 - elev - asin((Const.RE/(Const.RE + hI)).* cos(elev));
end

function isGnssSys = isGnssSys(sat_idx, idx_sys)
isGnssSys = sat_idx>=idx_sys(1) && sat_idx<=idx_sys(end);
end
