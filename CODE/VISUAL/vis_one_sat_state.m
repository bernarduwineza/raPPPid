% Plot the state of one (1 satellite used in the final solution
% This will use data saved in `storeData`, `model_save`, and `satellites`.

set(0,'DefaultFigureWindowStyle','docked')
gps_sats = [1,3, 7, 9, 11, 14, 18, 22, 23, 26, 31, 32];
gal_sats = 200+ [2, 7, 8, 13, 19, 21, 26, 27, 30, 33];


for sat = [1,3] %[gps_sats, gal_sats]
    plot_one_sat(storeData, model_save, settings, satellites, sat);
    axis tight;
end
plot_common_params(storeData, model_save, settings)

function plot_common_params(storeData, model_save, settings)
fig_title = sprintf('Common params for %s', settings.INPUT.file_obs(end-23:end));
togglefig(fig_title); clf;

clr = {'.b', 'r.', 'g.'};

num_rows = 6;
num_cols = 2;

t = tiledlayout(num_rows,num_cols,'TileSpacing','Compact');
axs = 1;
x_text = 'Elapsed time (sec)';
clf;

% epoch in elapsed GPS time
elapsed_time_sec = storeData.gpstime-storeData.gpstime(1);

% --- Plot Tropo corrections
tropo_corr   = model_save.zhd + model_save.zwd;
tropo_est    = storeData.zhd + storeData.zwd + storeData.param(:,4);
tropo_var    = sqrt(storeData.param_var(:,4));

ax(axs) = nexttile; axs = axs + 1;
title('Tropo Delay')
hold on
grid minor;
grid on;
xlabel(x_text);

ylabel('Tropo delay (m)');

plot(elapsed_time_sec, tropo_corr, clr{3});
plot(elapsed_time_sec, tropo_est, clr{1});
plot(elapsed_time_sec, tropo_est + tropo_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, tropo_est - tropo_var, 'r-.', 'LineWidth',2);
legend('Tropo correction', 'Tropo estimate')


% --- Plot receiver GPS (reference) clock error estimate and standard deviation
clk_err_gps     = storeData.param(:,5);
clk_err_gps(clk_err_gps==0)=NaN;
clk_err_var     = sqrt(storeData.param_var(:,5));

ax(axs) = nexttile; axs = axs + 1;
title('GPS receiver clock error')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Receiver clock error (m)');

plot(elapsed_time_sec, clk_err_gps, clr{1});
plot(elapsed_time_sec, clk_err_gps-clk_err_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, clk_err_gps+clk_err_var, 'r-.', 'LineWidth',2);

% --- Plot receiver GPS (reference) clock drift and standard deviation
clk_drift_gps     = storeData.clock_drift(:,1);
clk_drift_gps(clk_drift_gps==0)=NaN;
clk_drift_gps(1) = NaN;
clk_drift_var     = sqrt(storeData.clock_drift_var(:));

ax(axs) = nexttile; axs = axs + 1;
title('GPS receiver clock drift')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Receiver clock drift (m/s)');

plot(elapsed_time_sec, clk_drift_gps, clr{1});
plot(elapsed_time_sec, clk_drift_gps-clk_drift_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, clk_drift_gps+clk_drift_var, 'r-.', 'LineWidth',2);

% --- Plot receiver GAL clock estimate and standard deviation
clk_err_gal     = storeData.param(:,11);
clk_err_gal(clk_err_gal==0)=NaN;
clk_err_var     = sqrt(storeData.param_var(:,11));

ax(axs) = nexttile; axs = axs + 1;
title('GAL receiver clock offset')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Receiver clock error (m)');

plot(elapsed_time_sec, clk_err_gal, clr{1});
plot(elapsed_time_sec, clk_err_gal-clk_err_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, clk_err_gal+clk_err_var, 'r-.', 'LineWidth',2);

% --- Plot receiver DCB estimate
clk_dcb_gps     = storeData.param(:,6);
clk_dcb_gps(clk_dcb_gps==0)=NaN;
clk_dcb_var_gps = sqrt(storeData.param_var(:,6));

clk_dcb_gal     = storeData.param(:,12);
clk_dcb_gal(clk_dcb_gal==0)=NaN;
clk_dcb_var_gal = sqrt(storeData.param_var(:,12));

ax(axs) = nexttile; axs = axs + 1;
title('Receiver DCB for each system')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Receiver DCB (m)');

plot(elapsed_time_sec, clk_dcb_gps, clr{1});
plot(elapsed_time_sec, clk_dcb_gps - clk_dcb_var_gps, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, clk_dcb_gps + clk_dcb_var_gps, 'r-.', 'LineWidth',2);

plot(elapsed_time_sec, clk_dcb_gal, clr{1});
plot(elapsed_time_sec, clk_dcb_gal - clk_dcb_var_gal, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, clk_dcb_gal + clk_dcb_var_gal, 'r-.', 'LineWidth',2);

% legend('GPS', 'GAL')

% --- Link axes
linkaxes(ax, 'x')
end



function plot_one_sat(storeData, model_save, settings, satellites, sat)
fig_title = sprintf('Satellite %03d', sat);
togglefig(fig_title); clf;

% --- Some constants
lat_min = settings.IONO.Bspline.lat_min+1;
lat_max = settings.IONO.Bspline.lat_max;
lon_min = settings.IONO.Bspline.lon_min;
lon_max = settings.IONO.Bspline.lon_max;

H = 450e3;
K = settings.IONO.Bspline.K;
code_res_y_lim = [-10, 10];
phase_res_y_lim = [-.1, .1];
mp_y_lim = [-1, 1];
T_m  = 600;

clr = {'.b', 'r.', 'g.'};
% true if GNSS was processed and should be plotted
isGPS = settings.INPUT.use_GPS;
isGLO = settings.INPUT.use_GLO;
isGAL = settings.INPUT.use_GAL;
isBDS = settings.INPUT.use_BDS;

% epoch in elapsed GPS time
elapsed_time_sec = storeData.gpstime-storeData.gpstime(1);


% test satellite
test_sat = sat;

no_freqs = settings.INPUT.proc_freqs;

num_rows = 6;
num_cols = 2;

t = tiledlayout(num_rows,num_cols,'TileSpacing','Compact');
axs = 1;
x_text = 'Elapsed time (sec)';
clf;

%%% --- Plot satellite elevation
sat_el   = satellites.elev(:,test_sat);
sat_el(sat_el==0) = NaN;

ax(axs) = nexttile;
axs = axs + 1;
scatter(elapsed_time_sec, sat_el, [], sat_el, 's', 'filled');
title('Elevation')
grid minor;
grid on;
xlabel(x_text);
ylabel('Elevation (deg)');

elev_mask = 12;
len = 901;
len1 = elev_mask*10;
len2 = len-len1;
cbar1 = gray(len1*2);
cbar2 = jet(len2);
cm = [cbar1(end-1:-1:len1, :); cbar2];
colormap(cm)
clim([0,90])
box on;

%
% patch(elapsed_time_sec, sat_el, sat_el, 'EdgeColor','interp','Marker','.','MarkerFaceColor','flat');
colorbar;
clim([0 90]);

%%% --- Plot pseudorange time differences for all frequencies
if settings.INPUT.num_freqs > 1
    code_1   = storeData.C1(:,test_sat);
    code_1(code_1==0) = NaN;
    code_2   = storeData.C2(:,test_sat);
    code_2(code_2==0) = NaN;

    mp1 = storeData.mp1(:,test_sat);
    mp1(mp1==0) = NaN;

    mp2 = storeData.mp2(:,test_sat);
    mp2(mp2==0) = NaN;

    if(T_m~=0)
        mp1_tm = mp1 - movmean(mp1, T_m);
        mp2_tm = mp2 - movmean(mp2, T_m);
    else
        mp1_tm = mp1;
        mp2_tm = mp2;
    end 
else
    code_1   = model_save.code(:,test_sat);
    code_1(code_1==0) = NaN;
    code_2 = zeros(size(code_1));
end
% subplot(num_rows, num_cols, 1);
ax(axs) = nexttile;
axs = axs + 1;
% title('MP1')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('MP1 (m)');
box on;

% plot(code_1, '.')
% plot(code_2, '.')
if settings.INPUT.num_freqs > 1
    plot(elapsed_time_sec, mp1_tm, clr{1})
    ylim(mp_y_lim)
end


%%% --- Plot phase for all frequencies
if settings.INPUT.num_freqs > 1
    phase_1             = storeData.L1(:,test_sat);
    phase_1(phase_1==0) = NaN;
    phase_2             = storeData.L2(:,test_sat);
    phase_2(phase_2==0) = NaN;
else
    phase_1             = model_save.phase(:,test_sat);
    phase_1(phase_1==0) = NaN;
    phase_2 = zeros(size(phase_1));
end
% subplot(num_rows, num_cols, 2);
ax(axs) = nexttile; axs = axs + 1;
% title('MP2')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('MP2 (m)');
box on;

% plot(phase_1, '.')
% plot(phase_2, '.')
% plot(elapsed_time_sec, [NaN; NaN; diff(phase_1, 2)], clr{1})
% if settings.INPUT.num_freqs > 1
%     plot(elapsed_time_sec, [NaN; NaN; diff(phase_2, 2)], clr{2})
% end
if settings.INPUT.num_freqs > 1
    plot(elapsed_time_sec, mp2_tm, clr{2})
end
ylim(mp_y_lim)

% --- Plot code residuals
code_res_1      = storeData.residuals_code_1(:,test_sat);
code_res_1(code_res_1==0)=NaN;
code_res_var_1  = sqrt(storeData.residuals_code_var_1(:,test_sat));

% subplot(num_rows, num_cols, 3);
ax(axs) = nexttile; axs = axs + 1;
title('Code residual 1 + STD')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Code residual (m)');
box on;

plot(elapsed_time_sec, code_res_1, clr{1});
plot(elapsed_time_sec, code_res_1-code_res_var_1, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, code_res_1+code_res_var_1, 'r-.', 'LineWidth',2);
ylim(code_res_y_lim)

if settings.INPUT.num_freqs > 1
    code_res_2      = storeData.residuals_code_2(:,test_sat);
    code_res_2(code_res_2==0) = NaN;
    code_res_var_2  = sqrt(storeData.residuals_code_var_2(:,test_sat));

    ax(axs) = nexttile; axs = axs + 1;
    title('Code residual 2 + STD')
    hold on
    grid minor;
    grid on;
    xlabel(x_text);
    ylabel('Code residual (m)');
    box on;

    plot(elapsed_time_sec, code_res_2, clr{1});
    plot(elapsed_time_sec, code_res_2-code_res_var_2, 'r-.', 'LineWidth',2);
    plot(elapsed_time_sec, code_res_2+code_res_var_2, 'r-.', 'LineWidth',2);
end
% legend('Code 1', 'Code 2')

ylim(code_res_y_lim)


% --- Plot phase residual
phase_res_1  = storeData.residuals_phase_1(:,test_sat);
phase_res_1(phase_res_1==0)=NaN;
phase_res_var_1     = sqrt(storeData.residuals_phase_var_1(:,test_sat));


% subplot(num_rows, num_cols, 3);
ax(axs) = nexttile; axs = axs + 1;
title('Phase 1 residual + STD')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Phase 1 residual (m)');
ylim(phase_res_y_lim);
box on;

plot(elapsed_time_sec, phase_res_1, clr{1});
plot(elapsed_time_sec, phase_res_1-phase_res_var_1, 'r.', 'LineWidth',2);
plot(elapsed_time_sec, phase_res_1+phase_res_var_1, 'r.', 'LineWidth',2);

if settings.INPUT.num_freqs > 1
    phase_res_2  = storeData.residuals_phase_2(:,test_sat);
    phase_res_2(phase_res_2==0)=NaN;
    phase_res_var_2     = sqrt(storeData.residuals_phase_var_2(:,test_sat));

    ax(axs) = nexttile; axs = axs + 1;
    title('Phase 2 residual + STD')
    hold on
    grid minor;
    grid on;
    xlabel(x_text);
    ylabel('Phase 2 residual (m)');
    ylim(phase_res_y_lim);
    box on;

    plot(elapsed_time_sec, phase_res_2, clr{1});
    plot(elapsed_time_sec, phase_res_2-phase_res_var_2, 'r.', 'LineWidth',2);
    plot(elapsed_time_sec, phase_res_2+phase_res_var_2, 'r.', 'LineWidth',2);
    % legend('Phase 1', 'Phase 2')
end


% --- Plot Iono corrections
if settings.INPUT.num_freqs > 1
    iono_corr   = model_save.iono(:,test_sat,1);
else
    iono_corr   = model_save.iono(:,test_sat);
end
iono_corr(iono_corr==0)=NaN;
iono_est    = storeData.iono_est(:,test_sat);
iono_est(iono_est==0)=NaN;

% - compute iono var
iono_var = zeros(size(iono_est));
for i =1:length(storeData.bspline_coeff)
    idx_sats = storeData.lat_pp(i,:) ~= 0;
    num_params = length(storeData.param_sigma{i});
    phi_rx  = storeData.posFloat_geo(i,1);
    lam_rx  = storeData.posFloat_geo(i,2);
    Az      = satellites.az(i,test_sat);
    El      = satellites.elev(i,test_sat);
    [lat_pp, lon_pp] = calcIPP(phi_rx, lam_rx, Az*pi/180, El*pi/180, H);
    lat_lon_bases = compute_bspline_bases(settings, lat_pp, lon_pp);

    % idx_iono    = length(storeData.param(1,:))-settings.IONO.Bspline.K+1:length(storeData.param(1,:));
    d_iono      = storeData.bspline_coeff(i,:)';
    P_iono      = storeData.param_sigma{i};
    P_iono      = P_iono(num_params-K:num_params-1,num_params-K:num_params-1);
    iono_var(i) = sqrt(lat_lon_bases'*P_iono*lat_lon_bases);
end

ax(axs) = nexttile; axs = axs + 1;
title('Iono Delay (Correction & Estimate)')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Iono delay (m)');
box on; 

plot(elapsed_time_sec, iono_corr, clr{3});
plot(elapsed_time_sec, iono_est, clr{1});
plot(elapsed_time_sec, iono_est-iono_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, iono_est+iono_var, 'r-.', 'LineWidth',2);
legend('Iono correction', 'Iono estimate')


% --- Plot ambiguities
N1_amb     = storeData.N_1(:,test_sat);
N1_amb(N1_amb==0)=NaN;
N1_amb_var = sqrt(storeData.N_var_1(:,test_sat));


ax(axs) = nexttile; axs = axs + 1;
title('N_1 & N_2 float ambiguities')
hold on
grid minor;
grid on;
xlabel(x_text);
ylabel('Ambiguity (cy)');
box on;

plot(elapsed_time_sec, N1_amb+N1_amb_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, N1_amb-N1_amb_var, 'r-.', 'LineWidth',2);
plot(elapsed_time_sec, N1_amb, clr{1});

if settings.INPUT.num_freqs > 1
    N2_amb     = storeData.N_2(:,test_sat);
    N2_amb(N2_amb==0)=NaN;
    N2_amb_var = sqrt(storeData.N_var_2(:,test_sat));

    plot(elapsed_time_sec, N2_amb+N2_amb_var, 'r-.', 'LineWidth',2);
    plot(elapsed_time_sec, N2_amb-N2_amb_var, 'r-.', 'LineWidth',2);
    plot(elapsed_time_sec, N2_amb, clr{1});
end

% --- Link axes
linkaxes(ax, 'x')

end