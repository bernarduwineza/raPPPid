function [AgentData] = InitializeAgentProcessing(settings)

    %% -+-+-+-+-+-PREPARATIONS-+-+-+-+-+-
% Create path of results folder, the folder itself is created after the processing
if exist('settings', 'var')     % PPP_main.m started from GUI or with settings as input
    settings.PROC.output_dir = createResultsFolder(settings.INPUT.file_obs, settings.PROC.name);
else                            % PPP_main.m started without settings as input
    [AgentData.FileName, AgentData.PathName] = uigetfile('*.mat', 'Select a Settings File', Path.RESULTS');
    if ~FileName;       return;         end         % stop if no file selected
    AgentData.PathName = relativepath(AgentData.PathName);              % convert absolute path to relative path
    load([AgentData.PathName, '/', AgentData.FileName], 'settings');    % load selected settings.mat-file
    settings.PROC.output_dir = createResultsFolder(settings.INPUT.file_obs, settings.PROC.name);
end

AgentData.tStart = tic; 
warning off; 
fclose all;

% check if output should be printed to command window
AgentData.bool_print = ~settings.INPUT.bool_parfor;

% Prepare waitbar and print out of epochs to command window 
if AgentData.bool_print
    AgentData.WBAR = waitbar(0, 'Reading data, wait for processing start...', 'Name', 'Preparing...');
end

% Initialization and definition of global variable
global STOP_CALC;	% for stopping calculations before last epoch, set to zero in GUI_PPP

% number of processed frequencies
[settings.INPUT.num_freqs, settings.INPUT.proc_freqs] = CountProcessedFrequencies(settings);

% number of processed GNSS
settings.INPUT.use_GNSS = settings.INPUT.use_GPS + settings.INPUT.use_GLO + settings.INPUT.use_GAL + settings.INPUT.use_BDS;

% check which biases are applied
AgentData.bool_sinex = strcmp(settings.BIASES.phase(1:3), 'WHU') || ...
    any(strcmp(settings.BIASES.code, {'CAS Multi-GNSS DCBs','CAS Multi-GNSS OSBs','DLR Multi-GNSS DCBs','CODE OSBs','CNES OSBs','CODE MGEX','WUM MGEX','CNES MGEX','GFZ MGEX','CNES postprocessed'}));
AgentData.bool_manual_sinex = strcmp(settings.BIASES.code, 'manually') && settings.BIASES.code_manually_Sinex_bool;
AgentData.bool_CODE_dcb = strcmp(settings.BIASES.code, 'CODE DCBs (P1P2, P1C1, P2C2)') || ... 
    (strcmp(settings.BIASES.code, 'manually') && settings.BIASES.code_manually_DCBs_bool);
AgentData.bool_CNES_archive = settings.ORBCLK.bool_precise~=1   &&   strcmpi(settings.ORBCLK.CorrectionStream, 'CNES Archive');
AgentData.bool_brdc_TGD = strcmp(settings.BIASES.code, 'Broadcasted TGD');



%% -+-+-+-+-+-READ/PREPARE INPUT DATA-+-+-+-+-+-
AgentData.q_update = 15;      % [epochs], update rate of, for example, waitbar 
% Read Input Data from Files (Ephemerides, Station Data, etc.)
[AgentData.input, AgentData.obs, settings, AgentData.RINEX] = readAllInputFiles(settings);


% if necessary convert time frame to epochs of RINEX observation file 
% ||| missing data epochs are not considered
if settings.PROC.timeSpan_format_epochs
    settings.PROC.epochs(1) = floor(settings.PROC.timeFrame(1));    % round to be on the safe side
    settings.PROC.epochs(2) = ceil (settings.PROC.timeFrame(2));
elseif settings.PROC.timeSpan_format_SOD
    settings.PROC.epochs = sod2epochs(AgentData.RINEX, AgentData.obs.epochheader, settings.PROC.timeFrame, AgentData.obs.rinex_version);  
elseif settings.PROC.timeSpan_format_HOD
    settings.PROC.epochs = sod2epochs(AgentData.RINEX, AgentData.obs.epochheader, settings.PROC.timeFrame*3600, AgentData.obs.rinex_version);   % *3600 in order to convert from seconds to hours
end
% convert start of fixing from minutes (GUI) to epochs for calculations
if settings.AMBFIX.bool_AMBFIX
    settings.AMBFIX.start_WL = settings.AMBFIX.start_WL_sec/AgentData.obs.interval;
    settings.AMBFIX.start_NL = settings.AMBFIX.start_NL_sec/AgentData.obs.interval;
    settings.AMBFIX.start_fixing = [settings.AMBFIX.start_WL, settings.AMBFIX.start_NL];
end


    
%% -+-+-+-+-+-INITIALIZATION OF PROCESSING -+-+-+-+-+-

% create Epoch, satellites, storeData and update obs
[AgentData.Epoch, AgentData.satellites, AgentData.storeData, AgentData.obs, AgentData.model_save, AgentData.Adjust] = initProcessing(settings, AgentData.obs);

% create init_ambiguities
AgentData.init_ambiguities = NaN(3, 399);     % columns = satellites, rows = frequencies

% initialize matrices for Hatch-Melbourne-W�bbenab (HMW) LC between
% different frequencies, e.g., for GPS L1-L2-L5 processing
AgentData.HMW_12 = []; 
AgentData.HMW_23 = []; 
AgentData.HMW_13 = [];
if settings.AMBFIX.bool_AMBFIX
    AgentData.HMW_12 = zeros(settings.PROC.epochs(2)-settings.PROC.epochs(1)+1, 399);  % e.g., Wide-Lane (WL)
    AgentData.HMW_23 = zeros(settings.PROC.epochs(2)-settings.PROC.epochs(1)+1, 399);  % e.g., Extra-Wide-Lane (EW)
    AgentData.HMW_13 = zeros(settings.PROC.epochs(2)-settings.PROC.epochs(1)+1, 399);  % e.g., Medium-Lane (ML)
end

% Creating q for loop of epoch calculations, q is used for indexing the epochs
AgentData.q_range = 1:settings.PROC.epochs(2)-settings.PROC.epochs(1)+1; % one epoch more than the settings in GUI

% time which is needed for reading all input data and going to start-epoch
AgentData.read_time = toc(AgentData.tStart);

% check if P-code is processed for GPS and Glonass
% ||| does not work for P1 on 2nd freq or P2 on 1st freq
if AgentData.bool_CODE_dcb
    idx_c = AgentData.obs.use_column{1, 4};
    AgentData.bool_P1_GPS = strcmp('P1', AgentData.obs.types_gps(2*idx_c-1:2*idx_c));
    idx_c = AgentData.obs.use_column{1, 5};
    AgentData.bool_P2_GPS = strcmp('P2', AgentData.obs.types_gps(2*idx_c-1:2*idx_c));
    idx_c = AgentData.obs.use_column{2, 4};
    AgentData.bool_P1_GLO = strcmp('P1', AgentData.obs.types_glo(2*idx_c-1:2*idx_c));
    idx_c = AgentData.obs.use_column{2, 5};
    AgentData.bool_P2_GLO = strcmp('P2', AgentData.obs.types_glo(2*idx_c-1:2*idx_c));
end

% check for which satellites no precise orbits or clocks exist
if settings.ORBCLK.bool_precise
    settings = checkPreciseOrbitClock(settings, AgentData.input);
end

% Print stuff for epoch 0
if AgentData.bool_print
    AgentData.bspace = char(8,8,8,8,8,8,8,8,8,8,8,8,8,8);     % backspaces for printing actual epoch
    estr = sprintf('Epoch %d\n', 0);
    fprintf('%s',estr);
    AgentData.l_estr = length(estr);
    if ishandle(AgentData.WBAR)
        AgentData.WBAR.Name = ['Processing ' AgentData.obs.stationname sprintf(' %04.0f', AgentData.obs.startdate(1)) sprintf('/%03.0f', AgentData.obs.doy)];
    end
end

% Create agent data struct containing Epoch, satellites, storeData, obs, model_save, Adjust, settings
AgentData.settings          = settings;