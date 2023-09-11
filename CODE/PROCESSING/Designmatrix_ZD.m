function [Adjust, Epoch] = Designmatrix_ZD(Adjust, Epoch, model, settings)
% Create Designmatrix A and observed minus computed vector omc for code 
% and phase solution for Zero-Difference-Model
% 
% INPUT:
%   Adjust      ...
% 	Epoch       ...
% 	model           struct model
%   settings        settings from GUI
% OUTPUT: 
%   Adjust      updated with A and omc
%
% Revision:
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************


%% Preparations
% get variables from settings
num_freq = settings.INPUT.proc_freqs;   % number of processed frequencies
% get variables from Adjust
param   = Adjust.param;         % parameter estimation from last epoch
% get variables from Epoch
obs_code  = Epoch.code;         % code observations
obs_phase = Epoch.phase;        % phase observations
no_sats = length(Epoch.sats);  	% number of satellites of current epoch
isGPS = repmat(Epoch.gps, num_freq, 1); 	% logical vector for all frequencies, true if GPS satellite
isGLO = repmat(Epoch.glo, num_freq, 1);     % logical vector for all frequencies, true if GLO satellite
isGAL = repmat(Epoch.gal, num_freq, 1);  	% logical vector for all frequencies, true if GAL satellite
isBDS = repmat(Epoch.bds, num_freq, 1);  	% logical vector for all frequencies, true if BDS satellite


NO_PARAM = DEF.NO_PARAM_ZD;         % 8 estimated parameters
s_f = no_sats*num_freq;             % satellites x frequencies
inf_val = 1e8;                      % infinite value
cutoff = Epoch.exclude(:);           % satellite under cutoff angle?
usePhase = ~Epoch.cs_found(:);      % use phase observation if no CS was found
rho = model.rho(:);                 % geometric distance
sat_x = repmat(model.Rot_X(1,:)', 1, num_freq);   % satellite ECEF position x
sat_y = repmat(model.Rot_X(2,:)', 1, num_freq);   % satellite ECEF position y
sat_z = repmat(model.Rot_X(3,:)', 1, num_freq);   % satellite ECEF position z


%% Create Design Matrix and observed-minus-computed
% Design matrix [2 * #sats * #frequencies x #estimated parameters + #sats * #frequencies]
A   = zeros(2*s_f, NO_PARAM+s_f);
% observed minus computed observation row vector [2*number of sats*number of frequencies]
omc = zeros(2*s_f,1);                  
code_row = 1:2:2*s_f;   	% rows for code  obs [1,3,5,7,...]
phase_row = 2:2:2*s_f;  	% rows for phase obs [2,4,6,8,...]

% --- observed minus computed
omc(code_row,1)	 = (obs_code(:)  - model.model_code(:))  .*  ~cutoff; 	% for code-observations
omc(phase_row,1) = (obs_phase(:) - model.model_phase(:)) .*  ~cutoff .*  usePhase;    % for phase-observations
% usePhase = true, satellite is tracked long enough
% ~cutoff = true, satellite is true

% --- Partial derivatives
% coordinates
dR_dx    = -(sat_x(:)-param(1)) ./ rho; 	% x
dR_dy    = -(sat_y(:)-param(2)) ./ rho; 	% y
dR_dz    = -(sat_z(:)-param(3)) ./ rho; 	% z

% zenith wet delay
dR_dtrop = model.mfw(:); 	% get wet troposphere mapping-function

% initialize matrices for receiver clock error, time offsets and DCBs
dR_time_dcb = zeros(s_f, 12);
% receiver clock error
if settings.INPUT.use_GPS
    dR_time_dcb(:, 1) = isGPS + isGLO + isGAL + isBDS;       % gps
end
% time offsets
dR_time_dcb(:, 4) = isGLO;       % glonass
dR_time_dcb(:, 7) = isGAL;       % galileo
dR_time_dcb(:,10) = isBDS;       % beidou
% Differential Code Biases
if settings.BIASES.estimate_rec_dcbs
    % DCBs between 1st and 2nd frequency
    if num_freq > 1
        frq_2 = (no_sats+1):2*no_sats;
        dR_time_dcb(frq_2, 2) = -isGPS(frq_2);   	% gps
        dR_time_dcb(frq_2, 5) = -isGLO(frq_2);   	% glonass
        dR_time_dcb(frq_2, 8) = -isGAL(frq_2);   	% galileo
        dR_time_dcb(frq_2,11) = -isBDS(frq_2);   	% beidou
    end
    % DCBs between 1st and 3rd frequency
    if num_freq > 2
        frq_3 = (2*no_sats+1):3*no_sats;
        dR_time_dcb(frq_3, 3) = -isGPS(frq_3);  	% gps
        dR_time_dcb(frq_3, 6) = -isGLO(frq_3);   	% glonass
        dR_time_dcb(frq_3, 9) = -isGAL(frq_3);    	% galileo
        dR_time_dcb(frq_3,12) = -isBDS(frq_3);   	% beidou
    end
end

% Ambiguities
amb_c = zeros(s_f,s_f);     % ambiguity-part of A-Matrix for code
amb_p = eye(s_f,s_f);       % ambiguity N expressed in meters

% --- Build A-Matrix
A(code_row,:)  = [dR_dx, dR_dy, dR_dz, dR_dtrop, dR_time_dcb, amb_c] .*  ~cutoff;
A(phase_row,:) = [dR_dx, dR_dy, dR_dz, dR_dtrop, dR_time_dcb, amb_p] .*  ~cutoff .* usePhase;


%% add ionosphere estimation part to A and omc
if (strcmpi(settings.IONO.model,'Estimate with ... as constraint') ...
        || strcmpi(settings.IONO.model,'Estimate'))
    % --- Design Matrix A
    dR_diono_code_f1  =  Epoch.f1.^2 ./ Epoch.f1.^2;
    A_iono = diag(dR_diono_code_f1);
    if num_freq > 1      % 2nd frequency is processed
        dR_diono_code_f2  =  Epoch.f1.^2 ./ Epoch.f2.^2;
        A_iono_2 = diag(dR_diono_code_f2);
        A_iono = [A_iono; A_iono_2];
        if num_freq > 2      % 3rd frequency is processed
            dR_diono_code_f3  =  Epoch.f1.^2 ./ Epoch.f3.^2;
            A_iono_3 = diag(dR_diono_code_f3);
            A_iono = [A_iono; A_iono_3];
        end
    end
    A_iono = kron(A_iono,ones(2,1));                % duplicate for phase observation
    phase_rows = 2:2:size(A_iono,1);                % rows of phase observations
    A_iono(phase_rows,:) = -A_iono(phase_rows,:); 	% change sign for phase observations
    % Put Design-Matrix together
    A = [A, A_iono];
    
    if strcmpi(settings.IONO.model,'Estimate with ... as constraint') && Adjust.constraint
        % --- ionospheric Pseudo-observations
        A_iono_observ = [zeros(no_sats, NO_PARAM + s_f), eye(no_sats)];
        A = [A; A_iono_observ];
        % --- observed-minus-computed
        n = numel(Adjust.param);    % initialize estimated ionospheric delay
        iono_est = Adjust.param(n-no_sats+1:n);
        omc_iono = model.iono(:,1) - iono_est;
        omc = [omc(:); omc_iono(:)];
    end
end

%% add VTEC estimation part to A and omc
if (strcmpi(settings.IONO.model,'Estimate VTEC'))

    kj1= settings.IONO.Bspline.kj1;
    kj2= settings.IONO.Bspline.kj2;

    lat_min = settings.IONO.Bspline.lat_min;
    lat_max = settings.IONO.Bspline.lat_max;
    lon_min = settings.IONO.Bspline.lon_min;
    lon_max = settings.IONO.Bspline.lon_max;

    % Iono height
    H = 450e3;
    % Receiver position in geodetic coordinates
    pos_geo = cart2geo(Adjust.param(1:3));

    % Iono pierce point
    [lat_pp, lon_pp] = calcIPP(pos_geo.ph, pos_geo.la, model.az*pi/180, model.el*pi/180, H);
    
    % Basis functions for receiver pierce points
    [lat_lon_bases] = compute_bspline_bases(settings, lat_pp(:,1), lon_pp(:,1));

    % Control points from the GIM model
    lat_ctrl_pts = lat_min:5:lat_max;
    lon_ctrl_pts = lon_min:5:lon_max;

    % --- Ttr....transmission time/time of emission
    % code_dist = Epoch.code(i_sat);            % before 14.1.2021
    code_dist = mean(Epoch.code(1,:), 'omitnan');   % should be more stable
    tau = code_dist/Const.C;    % approximate signal runtime from sat. to rec.
    % time of emission [sow (seconds of week)] = time of obs. - runtime
    Ttr = Epoch.gps_time - tau;       	
   
    % interpolate VTEC
    [lon_grid,lat_grid] = ndgrid(lon_ctrl_pts, lat_ctrl_pts);
    ctrl_pts = [lat_grid(:),lon_grid(:)];
    vtec_ctrl_pt_gim = zeros(length(ctrl_pts),1);
    
    ionex = read_ionex_TUW(settings.IONO.file_ionex);
    for ii = 1:length(ctrl_pts)
        vtec_ctrl_pt_gim(ii) = iono_gims(ctrl_pts(ii,1), ctrl_pts(ii,2), Ttr, ionex, settings.IONO.interpol);
    end

    % Basis functions' values for the control points
    lat_lon_bases_ctrl_pts = compute_bspline_bases(settings, ctrl_pts(:,1)*pi/180, ctrl_pts(:,2)*pi/180);
    lat_lon_bases_ctrl_pts(isnan(lat_lon_bases_ctrl_pts)) = 0; 

    % --- save the basis functions for the current epoch 
    Epoch.lat_lon_bases = lat_lon_bases; 
    Epoch.lat_pp = lat_pp(:,1)*180/pi;
    Epoch.lon_pp = lon_pp(:,1)*180/pi; 
 
    % --- Design Matrix A
    A_iono = model.iono_mf .* 40.3e16./Epoch.f1.^2 .* lat_lon_bases';
    % A_iono = A_iono.*dR_diono_code_f1; %%%DEBUG

    if num_freq > 1      % 2nd frequency is processed
        A_iono_2 = model.iono_mf .* 40.3e16./Epoch.f2.^2 .* lat_lon_bases';
        % A_iono_2 = A_iono_2./dR_diono_code_f2; %%%DEBUG
        A_iono = [A_iono; A_iono_2];
        if num_freq > 2      % 3rd frequency is processed
            A_iono_3 = model.iono_mf .* 40.3e16./Epoch.f3.^2 .* lat_lon_bases';
            % A_iono_3 = A_iono_3./dR_diono_code_f3; %%%DEBUG
            A_iono = [A_iono; A_iono_3];
        end
    end
    A_iono = kron(A_iono,ones(2,1));                % duplicate for phase observation
    phase_rows = 2:2:size(A_iono,1);                % rows of phase observations
    A_iono(phase_rows,:) = -A_iono(phase_rows,:); 	% change sign for phase observations

    % Remove rows with no phase or code measurements
    A_iono(code_row,:) = A_iono(code_row,:) .* ~cutoff;
    A_iono(phase_row,:) = A_iono(phase_row,:) .* ~cutoff .* usePhase;

    % Put Design-Matrix together
    A = [A, A_iono];

    % --- ionospheric Pseudo-observations
    A_iono_observ = [zeros(no_sats, NO_PARAM + s_f), lat_lon_bases'];
    A = [A; A_iono_observ];

    % --- observed-minus-computed
    n = numel(Adjust.param);    % initialize estimated ionospheric delay
    iono_vtec_est = lat_lon_bases' * Adjust.param(n-(kj1*kj2)+1:n);
    %---DEBUG_CLK
    iono_vtec_est = lat_lon_bases' * Adjust.param(n-(kj1*kj2):n-1);
    %---
    missing_vtec = iono_vtec_est == 0; 
    iono_vtec_est(missing_vtec) = model.iono_vtec(missing_vtec); % initialize missing VTEC values
    omc_iono_vtec = model.iono_vtec - iono_vtec_est;
    omc = [omc(:); omc_iono_vtec(:)];
    
    use_ctrl_pt = 1; 
    if use_ctrl_pt
        % --- observed-minus-computed for the control points
        iono_vtec_ctrl_est = lat_lon_bases_ctrl_pts' * Adjust.param(n-(kj1*kj2)+1:n);
        %---DEBUG_CLK
        iono_vtec_ctrl_est = lat_lon_bases_ctrl_pts' * Adjust.param(n-(kj1*kj2):n-1);
        %---
        % missing_vtec = iono_vtec_ctrl_est == 0;
        % iono_vtec_ctrl_est(missing_vtec) = model.iono_vtec(missing_vtec); % initialize missing VTEC values
        omc_iono_ctrl_vtec = vtec_ctrl_pt_gim - iono_vtec_ctrl_est;
        omc = [omc(:); omc_iono_ctrl_vtec(:)];

        % --- ionospheric Pseudo-observations
        A_iono_observ_ctrl = [zeros(length(vtec_ctrl_pt_gim), NO_PARAM + s_f), lat_lon_bases_ctrl_pts'];
        A = [A; A_iono_observ_ctrl];

    end
end

%---DEBUG_CLK: Add clock drift state 
dR_rx_clock_drift = zeros(size(A(:,1)));
A  = [A, dR_rx_clock_drift];
%---

%% save in Adjust
Adjust.A = A;
Adjust.omc = omc;


