function [] = vis_iono_comparison_plot(settings, storeData, xaxis_label, hours, resets, bool_obs, rgb)
    % Plot of modelled correction of ionospheric delay, estimation of
    % ionospheric delay and histogramm of the difference between the modelled
    % and estimated ionospheric delay - all for the 1st frequency
    %
    % INPUT:
    %   settings        struct, settings of processing
    %   storeData       struct, data from processing
    %   xaxis_label 	string, label for x-axis
    %   hours           vector, time [h] since beginning of processing
    %   resets          vector, time [s] of resets
    %   bool_obs        boolean, true if satellites is observed in epoch
    %   rgb             colors for plotting
    % OUTPUT:
    %   []
    %
    % This function belongs to raPPPid, Copyright (c) 2023, J.-B. Uwineza
    % *************************************************************************
    
    
    %% Preparation
    bool_corr = strcmpi(settings.IONO.model,'Correct with ...')   ||   strcmpi(settings.IONO.model,'Estimate with ... as constraint');
    bool_est  = strcmpi(settings.IONO.model,'Estimate with ... as constraint') || ...
        strcmpi(settings.IONO.model,'Estimate') || ...
        strcmpi(settings.IONO.model,'Estimate VTEC');
    bool_fixed = settings.PLOT.fixed;
    
    % change default colors for plotting
    coleurs_default = get(groot,'defaultAxesColorOrder');           % save default colors for reset
    set(groot,'defaultAxesColorOrder',rgb)
    
    % no_rows = bool_corr + bool_est + (bool_est&&bool_corr);
    no_rows = 2;
    
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
    idx_bds = 301:size(bool_obs,2);
    
    % boolean matrix for each GNSS, true if satellite in this epoch was observed
    gps_obs = bool_obs(:,idx_gps);
    glo_obs = bool_obs(:,idx_glo);
    gal_obs = bool_obs(:,idx_gal);
    bds_obs = bool_obs(:,idx_bds);

    gamma_phi_gps = Const.GPS_F2^2 / (Const.GPS_F1^2 - Const.GPS_F2^2); 
    gamma_rho_gps = Const.GPS_F2^2 / (Const.GPS_F2^2 - Const.GPS_F1^2); 
    gamma_glo = Const.GLO_k2^2 / (Const.GLO_F1^2 - Const.GLO_k2^2);
    gamma_gal = Const.GAL_F5b^2 / (Const.GAL_F1^2 - Const.GAL_F5b^2);
    gamma_bds = Const.BDS_F2^2 / (Const.BDS_F1^2 - Const.BDS_F2^2);

    gps_prn = 26;
    gal_prn = 33; 
    
    fig_iono_comp = figure('Name','Ionospheric Delay Comparison','NumberTitle','off');
    i_plot = 1;
    % add customized datatip
    dcm = datacursormode(fig_iono_comp);
    datacursormode on
    set(dcm, 'updatefcn', @vis_customdatatip_h)    
    
    %% Plot: Estimation of Ionospheric Delay
    if bool_est
        if ~bool_fixed    	% estimation from float solution
            sol_str = 'Float';
            iono_est = full(storeData.iono_est);       
        else                % estimation from fixed solution
            sol_str = 'Fixed';
            iono_est = full(storeData.iono_fixed);      
        end
        if isGPS
            prn = gps_prn;
            iono_est_G = iono_est(:,idx_gps);
            iono_corr_G = storeData.iono_corr(:, idx_gps);
            iono_meas_phase_G = gamma_phi_gps * (storeData.L1(:,idx_gps) - storeData.L2(:,idx_gps));
            iono_meas_rho_G = gamma_rho_gps * (storeData.C2(:,idx_gps) - storeData.C1(:,idx_gps));

            max_iono_corr = max(iono_corr_G(:)); 
            subplot(no_rows, noGNSS, i_plot);
            i_plot=i_plot+1;
            plotIonoComparisons(iono_est_G, iono_corr_G, iono_meas_phase_G, iono_meas_rho_G, ...
                hours, resets, gps_obs, 'GPS', xaxis_label, sol_str, prn)
        end
        if isGLO
            iono_est_R = iono_est(:,idx_glo);
            iono_corr_R = storeData.iono_corr(:, idx_glo);
            iono_meas_phase_R = gamma_glo * (storeData.L1(:,idx_glo) - storeData.L2(:,idx_glo));
            iono_meas_rho_R = gamma_glo * (storeData.C2(:,idx_glo) - storeData.C1(:,idx_glo));

            max_iono_corr = max(iono_corr_R(:)); 
            subplot(no_rows, noGNSS, i_plot);     
            i_plot=i_plot+1;
         plotIonoComparisons(iono_est_R, iono_corr_R, iono_meas_phase_R, iono_meas_rho_R, ...
          hours, resets, glo_obs, 'Glonass', xaxis_label, sol_str, prn)
        end
        if isGAL
            prn = gal_prn;
            iono_est_E = iono_est(:,idx_gal);
            iono_corr_E = storeData.iono_corr(:, idx_gal);
            iono_meas_phase_E = gamma_gal * (storeData.L1(:,idx_gal) - storeData.L2(:,idx_gal));
            iono_meas_rho_E = gamma_gal * (storeData.C2(:,idx_gal) - storeData.C1(:,idx_gal));

            max_iono_corr = max(iono_corr_E(:)); 
            subplot(no_rows, noGNSS, i_plot);  
            i_plot=i_plot+1;
         plotIonoComparisons(iono_est_E, iono_corr_E, iono_meas_phase_E, iono_meas_rho_E, ...
          hours, resets, gal_obs, 'Galileo', xaxis_label, sol_str, prn)
        end
        if isBDS
            iono_est_C = iono_est(:,idx_bds);
            iono_corr_C = storeData.iono_corr(:, idx_bds);
            iono_meas_phase_C = gamma_bds * (storeData.L1(:,idx_bds) - storeData.L2(:,idx_bds));
            iono_meas_rho_C = gamma_bds * (storeData.C2(:,idx_bds) - storeData.C1(:,idx_bds));

            max_iono_corr = max(iono_corr_C(:)); 
            subplot(no_rows, noGNSS, i_plot);       
            i_plot=i_plot+1;
            plotIonoComparisons(iono_est_C, iono_corr_C, iono_meas_phase_C, iono_meas_rho_C, ...
          hours, resets, bds_obs, 'BeiDou', xaxis_label, sol_str, prn)
        end
    end

    %% Plot Differences of estimates across the two receivers
    if isGPS
        subplot(no_rows, noGNSS, i_plot);
        i_plot=i_plot+1;
        plotIonoDifferences(iono_corr_G, iono_meas_phase_G, iono_meas_rho_G, ...
            hours, 'GPS', xaxis_label, sol_str, gps_prn)
    end

    if isGAL
        subplot(no_rows, noGNSS, i_plot);
        i_plot=i_plot+1;
        plotIonoDifferences(iono_corr_E, iono_meas_phase_E, iono_meas_rho_E, ...
            hours, 'Galileo', xaxis_label, sol_str, gal_prn)
    end

    
    %% reset to default colors
    set(groot,'defaultAxesColorOrder',coleurs_default)
    
end
    
%% Auxiliary Functions
% For plotting the estimated ionospheric delay
function [] = plotIonoComparisons(iono_est, iono_corr, iono_meas_phi, iono_meas_rho, ...
     hours, resets, gnss_obs, gnss, xaxis_label, sol_str, prn)
    % select PRN to plot 
    % prn = 26;
    iono_meas_phi(iono_meas_phi(:,prn) == 0, prn) = NaN;
    iono_corr(iono_corr(:,prn) == 0, prn) = NaN;

    % Plot
    iono_est(iono_est==0) = NaN;
    hold on
    % plot(hours, iono_est(:,prn), 'r.'); 
    plot(hours, iono_corr(:,prn), 'b.');
    % plot(hours, iono_meas_phi(:,prn) + (iono_est(1,prn) - iono_meas_phi(1,prn)) , 'c.');
    plot(hours, iono_meas_phi(:,prn), 'c.');

    plot(hours, iono_meas_rho(:,prn), 'm.');
    % if ~isempty(resets); vline(resets, 'k:'); end	% plot vertical lines for resets
    % create legend which datatooltip needs
    sys = gnss2char(gnss);
    prn = strcat(sys, num2str(mod(prn',100), '%02.0f'));      % for legend
    legend('$I_{GIM}$', '$I_{\phi,1}$,  $I_{\phi,2}$', 'Interpreter', 'latex', 'fontsize', 14);
    % style
    % xlim([floor(hours(1)) hours(end)])
    % ylim([-20 20])
    title([sol_str ' iono delay estimation for ' prn])
    ylabel('Ionospheric Delay [m]')
    Grid_Xon_Yon();
    xlabel(xaxis_label)
end

function [] = plotIonoDifferences(iono_corr, iono_meas_phi, iono_meas_rho, ...
     hours, gnss, xaxis_label, sol_str, prn)
    % First indices on all available reciever data
    first_idx = find(hours == 0); 
    iono_meas_phi(iono_meas_phi(:,prn)==0, prn) = NaN;
    iono_corr(iono_corr(:,prn) == 0, prn) = NaN;
    
    % shortest duration of data for each receiver
    if (length(first_idx) == 1)
        min_duration = length(iono_meas_phi)-1;
    else
        min_duration = min([length(iono_meas_phi) - max(first_idx), max(first_idx)]);
    end 

    % isolate iono estimates from all receivers
    iono_rec_phi_est = sparse(length(first_idx), min_duration + 1);
    iono_rec_rho_est = sparse(length(first_idx), min_duration + 1);
    for i=1:length(first_idx)
        iono_rec_phi_est(i, :) = iono_meas_phi(first_idx(i): first_idx(i) + min_duration,prn);
        iono_rec_rho_est(i, :) = iono_meas_rho(first_idx(i): first_idx(i) + min_duration,prn);
    end

    % average of receiver data for ~ 4hours
    avg_iono_meas_phi = mean(iono_rec_phi_est, 1);

    % plot the difference of the iono estimated from phase observations 
    hold on
    plot(hours(1: min_duration + 1), diff(iono_rec_phi_est, 1, 1), 'r.')
    plot(hours(1: min_duration + 1), -diff(iono_rec_phi_est, 1, 1) + (diff(iono_rec_phi_est(:,1), 1, 1)), 'm.')
    plot(hours(1: min_duration + 1), iono_corr(1:min_duration + 1, prn) - avg_iono_meas_phi' - ((iono_corr(1,prn) - avg_iono_meas_phi(1,1))) ,'b.')

    sys = gnss2char(gnss);
    prn = strcat(sys, num2str(mod(prn',100), '%02.0f'));      % for legend 
    legend('$I_{\phi,1} - I_{\phi,2}$', '$I_{\phi,1} - I_{\phi,2} \,- $ offset', '$\frac{1}{2}(I_{\phi,1} + I_{\phi,2})$ - $I_{GIM} \, -$ offset', ...
        'Interpreter', 'latex', 'fontsize', 14);
    % style
    % xlim([floor(hours(1)) hours(end)])
    % ylim([-20 20])
    title(['Iono delay differences for ' prn])
    ylabel('Ionospheric Delay Differences [m]')
    Grid_Xon_Yon();
    xlabel(xaxis_label)
end 