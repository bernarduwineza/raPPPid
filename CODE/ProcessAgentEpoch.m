function [mainAgentData] = ProcessAgentEpoch(mainAgentData, q)

% save data from last epoch (step q-1) and reset mainAgentData.Epoch
mainAgentData.Epoch.old = mainAgentData.Epoch;
mainAgentData.Epoch.q = q;            % save number of epoch
[mainAgentData.Epoch] = EpochlyReset_Epoch(mainAgentData.Epoch);

n = mainAgentData.settings.PROC.epochs(1) + q - 1;
if ~mainAgentData.settings.INPUT.bool_realtime
    % post-processing: check if end of file is reached
    if n > length(mainAgentData.obs.epochheader)
        mainAgentData.settings.PROC.epochs(2) = mainAgentData.settings.PROC.epochs(1) + q - 2;
        q = q - 1;          %#ok<FXSET>
        mainAgentData.Epoch.q = q;
        errordlg({['Epoch: ' sprintf(' %.0f',q)], 'End of Observation File reached!'}, ...
            [mainAgentData.obs.stationname sprintf(' %04.0f', mainAgentData.obs.startdate(1)) sprintf('/%03.0f',mainAgentData.obs.doy)]);
        return;
    end
end

% ----- read in epoch data from RINEX observation file -----
[mainAgentData.Epoch] = RINEX2epochData(mainAgentData.RINEX, mainAgentData.obs.epochheader, mainAgentData.Epoch, n, mainAgentData.obs.no_obs_types, mainAgentData.obs.rinex_version, mainAgentData.settings);
if ~mainAgentData.Epoch.usable
    if mainAgentData.bool_print
        fprintf('Epoch %d is skipped (not usable based on RINEX header)            \n', q)
    end
    [mainAgentData.Epoch, mainAgentData.storeData, mainAgentData.Adjust] = SkipEpoch(mainAgentData.Epoch, mainAgentData.storeData, mainAgentData.Adjust);
    return;
end

% reset solution
[mainAgentData.Adjust, mainAgentData.Epoch, mainAgentData.settings, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13, mainAgentData.storeData, mainAgentData.init_ambiguities] = ...
    resetSolution(mainAgentData.Adjust, mainAgentData.Epoch, mainAgentData.settings, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13, mainAgentData.storeData, mainAgentData.obs.interval, mainAgentData.init_ambiguities);

% ----- Find column of broadcast-ephemeris of current satellite -----
% check if satellites have broadcast ephemerides and are healthy, otherwise satellite is excluded
if mainAgentData.settings.ORBCLK.bool_brdc
    mainAgentData.Epoch = findEphCorr2Brdc(mainAgentData.Epoch, mainAgentData.input, mainAgentData.settings);
end

% ----- prepare observations -----
[mainAgentData.Epoch, mainAgentData.obs] = prepareObservations(mainAgentData.settings, mainAgentData.obs, mainAgentData.Epoch, q);

% ----- check, if enough satellites -----
bool_enough_sats = check_min_sats(mainAgentData.settings.INPUT.use_GPS, mainAgentData.settings.INPUT.use_GLO, mainAgentData.settings.INPUT.use_GAL, mainAgentData.settings.INPUT.use_BDS, ...
    sum(mainAgentData.Epoch.gps), sum(mainAgentData.Epoch.glo), sum(mainAgentData.Epoch.gal), sum(mainAgentData.Epoch.bds), mainAgentData.settings.INPUT.use_GNSS);
if ~bool_enough_sats
    if mainAgentData.bool_print
        fprintf('Less than %d usable satellites in epoch %d (%s)        \n', DEF.MIN_SATS, q, mainAgentData.Epoch.rinex_header);
    end
    [mainAgentData.Epoch, mainAgentData.storeData, mainAgentData.Adjust] = SkipEpoch(mainAgentData.Epoch, mainAgentData.storeData, mainAgentData.Adjust);
    return;
end

% ----- check, if epoch is excluded from processing -----
if ~isempty(mainAgentData.settings.PROC.excl_eps) && any(q == mainAgentData.settings.PROC.excl_eps)
    [mainAgentData.settings, mainAgentData.Epoch, mainAgentData.Adjust, mainAgentData.storeData, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13] = ...
        ExcludeEpoch(mainAgentData.settings, mainAgentData.Epoch, mainAgentData.Adjust, mainAgentData.storeData, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13, mainAgentData.bool_print);
    return;
end

% number of satellites in current epoch
mainAgentData.Epoch.no_sats = numel(mainAgentData.Epoch.sats);
% frequency
f1 = mainAgentData.Epoch.gps .* Const.GPS_F(strcmpi(DEF.freq_GPS_names, mainAgentData.settings.INPUT.gps_freq{1})) + mainAgentData.Epoch.gal .* Const.GAL_F(strcmpi(DEF.freq_GAL_names, mainAgentData.settings.INPUT.gal_freq{1})) + mainAgentData.Epoch.bds .* Const.BDS_F(strcmpi(DEF.freq_BDS_names, mainAgentData.settings.INPUT.bds_freq{1}));
f2 = mainAgentData.Epoch.gps .* Const.GPS_F(strcmpi(DEF.freq_GPS_names, mainAgentData.settings.INPUT.gps_freq{2})) + mainAgentData.Epoch.gal .* Const.GAL_F(strcmpi(DEF.freq_GAL_names, mainAgentData.settings.INPUT.gal_freq{2})) + mainAgentData.Epoch.bds .* Const.BDS_F(strcmpi(DEF.freq_BDS_names, mainAgentData.settings.INPUT.bds_freq{2}));
f3 = mainAgentData.Epoch.gps .* Const.GPS_F(strcmpi(DEF.freq_GPS_names, mainAgentData.settings.INPUT.gps_freq{3})) + mainAgentData.Epoch.gal .* Const.GAL_F(strcmpi(DEF.freq_GAL_names, mainAgentData.settings.INPUT.gal_freq{3})) + mainAgentData.Epoch.bds .* Const.BDS_F(strcmpi(DEF.freq_BDS_names, mainAgentData.settings.INPUT.bds_freq{3}));
f1(mainAgentData.Epoch.glo) = mainAgentData.Epoch.f1_glo;
f2(mainAgentData.Epoch.glo) = mainAgentData.Epoch.f2_glo;
f3(mainAgentData.Epoch.glo) = mainAgentData.Epoch.f3_glo;
mainAgentData.Epoch.f1 = f1;   mainAgentData.Epoch.f2 = f2;   mainAgentData.Epoch.f3 = f3;
% wavelength
lam1 = Const.C ./ f1;
lam2 = Const.C ./ f2;
lam3 = Const.C ./ f3;
mainAgentData.Epoch.l1 = lam1;
mainAgentData.Epoch.l2 = lam2;
mainAgentData.Epoch.l3 = lam3;
% get prn: GPS [0-99], Glonass [100-199], Galileo [200-299], BeiDou [300-399]
prn_Id = mainAgentData.Epoch.sats;
% increase epoch counter
mainAgentData.Epoch.tracked(prn_Id) = mainAgentData.Epoch.tracked(prn_Id) + 1;

% --- check SNR and signal strength threshold  ---
[mainAgentData.Epoch] = check_SNR(mainAgentData.Epoch, mainAgentData.settings, mainAgentData.obs.use_column);

% --- perform multipath detection ---
if mainAgentData.settings.OTHER.mp_detection
    [mainAgentData.Epoch] = checkMultipath(mainAgentData.Epoch, mainAgentData.settings, mainAgentData.obs.use_column, mainAgentData.obs.interval, mainAgentData.Adjust.reset_time);
end

% --- insert artificial cycle slip or multipath
% mainAgentData.Epoch = cycleSlip_articifial(mainAgentData.Epoch, obs.use_column);
% mainAgentData.Epoch = multipath_articifial(mainAgentData.Epoch, obs.use_column);

% --- perform Cycle-Slip detection ---
if contains(mainAgentData.settings.PROC.method,'Phase')   &&   (mainAgentData.settings.OTHER.CS.l1c1 || mainAgentData.settings.OTHER.CS.DF || mainAgentData.settings.OTHER.CS.Doppler || mainAgentData.settings.OTHER.CS.TimeDifference || mainAgentData.settings.PROC.LLI)
    [mainAgentData.Epoch] = cycleSlip(mainAgentData.Epoch, mainAgentData.settings, mainAgentData.obs.use_column);
end

% --- Adjust phase data to Code to limit the ambiguities ---
if strcmpi(mainAgentData.settings.PROC.method, 'Code + Phase') && mainAgentData.settings.PROC.AdjustPhase2Code
    [mainAgentData.init_ambiguities, mainAgentData.Epoch] = AdjustPhase2Code(mainAgentData.Epoch, mainAgentData.init_ambiguities, mainAgentData.Epoch.old.sats);
end

% --- get and apply satellite biases for current epoch ---
% get from correction stream (real-time code and phase biases)
if strcmp(mainAgentData.settings.ORBCLK.CorrectionStream,'manually') && (mainAgentData.settings.BIASES.code_corr2brdc_bool || mainAgentData.settings.BIASES.phase_corr2brdc_bool)
    mainAgentData.Epoch = apply_corr2brdc_biases(mainAgentData.Epoch, mainAgentData.settings, mainAgentData.input, mainAgentData.obs);
    % get from *.bia-file
elseif mainAgentData.bool_sinex || mainAgentData.bool_manual_sinex || mainAgentData.bool_CNES_archive
    mainAgentData.Epoch = apply_biases(mainAgentData.Epoch, mainAgentData.obs, mainAgentData.settings);
    % get from *.DCB-file
elseif mainAgentData.bool_CODE_dcb
    mainAgentData.Epoch = apply_DCBs(input, mainAgentData.settings, mainAgentData.Epoch, mainAgentData.bool_P1_GPS, mainAgentData.bool_P2_GPS, mainAgentData.bool_P1_GLO, mainAgentData.bool_P2_GLO);
    % get Broadcasted Time Group Delays
elseif mainAgentData.bool_brdc_TGD
    mainAgentData.Epoch = apply_TGDs(input, mainAgentData.settings, mainAgentData.Epoch);
end
% apply satellite biases [m]
mainAgentData.Epoch.C1 = mainAgentData.Epoch.C1 + mainAgentData.Epoch.C1_bias;
mainAgentData.Epoch.C2 = mainAgentData.Epoch.C2 + mainAgentData.Epoch.C2_bias;
mainAgentData.Epoch.C3 = mainAgentData.Epoch.C3 + mainAgentData.Epoch.C3_bias;
mainAgentData.Epoch.L1 = mainAgentData.Epoch.L1 + mainAgentData.Epoch.L1_bias;
mainAgentData.Epoch.L2 = mainAgentData.Epoch.L2 + mainAgentData.Epoch.L2_bias;
mainAgentData.Epoch.L3 = mainAgentData.Epoch.L3 + mainAgentData.Epoch.L3_bias;

% --- correct receiver biases ---
if mainAgentData.settings.INPUT.proc_freqs > 1 && ~mainAgentData.settings.BIASES.estimate_rec_dcbs
    mainAgentData.Epoch = correct_rec_biases(mainAgentData.Epoch, mainAgentData.obs);
end

% --- Build LCs and processed observations -> mainAgentData.Epoch.code/.phase ---
[mainAgentData.Epoch, mainAgentData.storeData] = create_LC_observations(mainAgentData.Epoch, mainAgentData.settings, mainAgentData.storeData, q);


% -+-+-+-+-+-START CALCULATION EPOCH-WISE SOLUTION-+-+-+-+-+-
[mainAgentData.Adjust, mainAgentData.Epoch, mainAgentData.model, mainAgentData.obs, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13] = ...
    ZD_processing(mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13, mainAgentData.Adjust, mainAgentData.Epoch, mainAgentData.settings, mainAgentData.input, mainAgentData.satellites, mainAgentData.obs);

% Save results from epoch-wise processing
[mainAgentData.satellites, mainAgentData.storeData, mainAgentData.model_save] = ...
    saveData(mainAgentData.Epoch, q, mainAgentData.satellites, mainAgentData.storeData, mainAgentData.settings, mainAgentData.Adjust, mainAgentData.model, mainAgentData.model_save, mainAgentData.HMW_12, mainAgentData.HMW_23, mainAgentData.HMW_13);
mainAgentData.Epoch.delta_windup = mainAgentData.model.delta_windup;        % [cycles], for calculation of Wind-Up correction in next epoch
