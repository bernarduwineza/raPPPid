function handles = StartCollaborativeProcessingFromGUI(handles)
% This function starts the Collaborative PPP processing from the GUI
%
% INPUT:
%	handles     struct, from raPPPid GUI
% OUTPUT:
%   results in RESULTS folder
%	handles     struct, from raPPPid GUI
%
% Revision:
%   ...
%
% This function belongs to raPPPid, Copyright (c) 2023, J.-B. Uwineza
% *************************************************************************

warning('off');

global STOP_CALC
STOP_CALC = 0;
l_start = ['Start: ' char(datetime('now'))];
bool_BATCH_PROC = handles.checkbox_batch_proc.Value;    % batch-processing enabled?
bool_parfor = handles.checkbox_parfor.Value;            % parfor loop enabled?


if bool_BATCH_PROC          % batch processing
    
    % reset plotting panel
    set(handles.edit_plot_path,'String', '');
    handles.paths.plotfile = '';
    un_check_plot_checkboxes(handles, 0);
    % reset true position into textfields of single plot
    set(handles.edit_x_true,  'String', '');
    set(handles.edit_y_true,  'String', '');
    set(handles.edit_z_true,  'String', '');
    % get data of batch processing table [cell]
    TABLE = GetTableData(handles.uitable_batch_proc.Data, 2, [6,9,12,15], [1 2], []);
    n = size(TABLE,1);        	% number of rows
    success = false(1,n);       % to save wich lines are successfully processed
    % get processing settings from GUI
    settings = getSettingsFromGUI(handles);
    % check if processing settings are valid in general
    valid_settings = checkProcessingSettings(settings, true);
    if ~valid_settings
        return;
    end
    if isempty(TABLE)
        errordlg({'Batch-Processing could not be started:'; 'Batch-Processing table is empty.'}, 'Fail')
        return
    end
    
    collabData = cell(n,1);
    sprintf('Collecting data for %d ...\n', n);
    t0 = tic;
    mainAgentData = cell(n,1);
    for i=1:n
        ROW_AGENT = TABLE(i,:);
        settings_agent = BatchProcessingPreparation(settings, ROW_AGENT);
        
        if checkProcessingSettings(settings_agent, false)
            collabData{i} = InitializeAgentProcessing(settings_agent);
            collabData{i}.agentNo = i;
        else
            sprintf('Invalid settings for %s\n', settings_agent.INPUT.file_obs);
            STOP_CALC = 0;
        end
        disp(i);
        if i==2
            mainAgentData{i} = collabData{i};
        end
    end
    mainAgentData =  mainAgentData(~cellfun('isempty', mainAgentData));
    mainAgentData = mainAgentData{1};
    mainAgentData.bool_print = 1;
    toc(t0)

    if mainAgentData.bool_print
        bspace = char(8,8,8,8,8,8,8,8,8,8,8,8,8,8);     % backspaces for printing actual epoch
        estr = sprintf('Epoch %d\n', 0);
        fprintf('%s',estr);
        l_estr = length(estr);
        % Prepare waitbar and print out of epochs to command window
        mainAgentData.WBAR = waitbar(0, 'Reading data, wait for processing start...', 'Name', 'Preparing...');
        mainAgentData.WBAR.Name = ['Processing ' mainAgentData.obs.stationname sprintf(' %04.0f', ...
            mainAgentData.obs.startdate(1)) sprintf('/%03.0f',mainAgentData.obs.doy)];
    end
    
    [~,file,ext] = fileparts(mainAgentData.settings.INPUT.file_obs);
    fprintf('\n---------------------------------------------------------------------\n');
    fprintf('%s%03d%s%03d%s\n%s\n\n','Batch Processing File #', 1, ' of ', n, ': ', [file ext])
    fprintf('\n---------------------------------------------------------------------\n');
    
    
    %% -+-+-+-+-+-START EPOCH PROCESSING-+-+-+-+-+-
    for q = mainAgentData.q_range
        
        if mainAgentData.bool_print       % printing out the number of the epoch
            estr = sprintf('\nEpoch %d\n', q);
            fprintf('%s%s', bspace(1:l_estr), estr);
            l_estr = length(estr);
        end
        
        % ----- COMPUTE SOLUTION FOR AN AGENT --------
        t0=tic;
        for i=1:n
            agentData = collabData{i};
            collabData{i} = ProcessAgentEpoch(agentData, q);
        end
        toc(t0)
        
        mainAgentData = collabData{2};
        % update waitbar
        if mainAgentData.bool_print && mod(q,mainAgentData.q_update) == 0 && ishandle(mainAgentData.WBAR)
            progress = q/mainAgentData.q_range(end);
            remains = (toc(mainAgentData.tStart)-mainAgentData.read_time)/q*mainAgentData.q_range(end) - ...
                (toc(mainAgentData.tStart)-mainAgentData.read_time); 	% estimated remaining time for epoch-wise processing in seconds
            mess = sprintf('%s%s\n%02.0f%s%d%s%d',...
                'Estimated remaining time: ', char(datetime(remains/(24*60*60), 'ConvertFrom', 'epochtime', 'Format', 'HH:mm:ss')), progress*100, '% processed, current epoch: ', q, ' of ', mainAgentData.q_range(end));
            waitbar(progress, mainAgentData.WBAR, mess)
        end
        
        % Stop processing: if user stopped manually or end of real-time processing
        if ~(mainAgentData.settings.INPUT.bool_batch && mainAgentData.settings.INPUT.bool_parfor)
            if mod(q,mainAgentData.q_update) == 0
                drawnow;
            end
            
            if STOP_CALC
                mainAgentData.settings.PROC.epochs(2) = q + (mainAgentData.settings.PROC.epochs(1)-1);
                break;
            end
        end
    end             % end of loop over epochs / epoch-wise calculations

    %% -+-+-+-+-+-VISUALS AND OUTPUT-+-+-+-+-+-
    
    mess = sprintf('Processing finished, opening plots...');
    if mainAgentData.bool_print
        % update waitbar
        if ishandle(mainAgentData.WBAR) 
            waitbar(1,mainAgentData.WBAR, mess);  
        end
    end
    
    % save epochs of reset
    mainAgentData.storeData.float_reset_epochs = mainAgentData.Adjust.float_reset_epochs;
    mainAgentData.storeData.fixed_reset_epochs = mainAgentData.Adjust.fixed_reset_epochs;
    
    % if processing was finished before the timespan defined in the GUI -> shrink variables
    [mainAgentData.satellites, mainAgentData.storeData, mainAgentData.model_save, mainAgentData.obs] = ...
        shrinkVariables(mainAgentData.satellites, mainAgentData.storeData, mainAgentData.model_save, mainAgentData.obs, mainAgentData.settings, q);
    
    % sparse variables
    [mainAgentData.satellites, mainAgentData.storeData, mainAgentData.model_save, mainAgentData.obs] = ...
        sparseVariables(mainAgentData.satellites, mainAgentData.storeData, mainAgentData.model_save, mainAgentData.obs, mainAgentData.settings);
    
    
    % create results folder
    mkdir(mainAgentData.settings.PROC.output_dir);
    
    % write settings from GUI into settings_summary.txt and settings.mat
    if mainAgentData.settings.EXP.settings_summary
        settings2txt(mainAgentData.settings, mainAgentData.obs, mainAgentData.input.OTHER.PCO.rec_error, mainAgentData.input.OTHER.PCV.rec_error, mainAgentData.tStart);
    end
    if mainAgentData.settings.EXP.settings
        settings = mainAgentData.settings;
        save(fullfile(mainAgentData.settings.PROC.output_dir, 'settings.mat'), '-struct', 'mainAgentData', 'settings')
    end
    
    % Create Output Files
    [mainAgentData.storeData] = create_output(mainAgentData.storeData, mainAgentData.obs, mainAgentData.settings, mainAgentData.Epoch.q );
    
    % write data4plot.mat
    if mainAgentData.settings.EXP.data4plot
        save([mainAgentData.settings.PROC.output_dir, '/data4plot.mat'], '-struct', 'mainAgentData', 'obs', 'settings');
        if mainAgentData.settings.EXP.storeData
            save([mainAgentData.settings.PROC.output_dir, '/data4plot.mat'], '-struct', 'mainAgentData', 'storeData', '-append')
        end
        if mainAgentData.settings.EXP.satellites
            save([settings.PROC.output_dir, '/data4plot.mat'], '-struct', 'mainAgentData', 'satellites', '-append')
        end
        if mainAgentData.settings.EXP.model_save
            save([settings.PROC.output_dir, '/data4plot.mat'], '-struct', 'mainAgentData', 'model_save', '-append')
        end
    end
    
    % push to workspace
    assignin('base',     'obs',         mainAgentData.obs        );
    assignin('base',     'settings',    mainAgentData.settings   );
    if mainAgentData.settings.EXP.satellites
        assignin('base', 'satellites',  mainAgentData.satellites );
    end
    if mainAgentData.settings.EXP.storeData
        assignin('base', 'storeData',   mainAgentData.storeData  );
        assignin('base', 'collabData',  collabData);
    end
    if settings.EXP.model_save
        assignin('base', 'model_save',  mainAgentData.model_save );
    end
    
    % close waitbar
    if mainAgentData.bool_print && ishandle(mainAgentData.WBAR)
        close(mainAgentData.WBAR);
    end
    
    % print current time
    fprintf('%s', char(datetime('now')))
    
    % Print final processing time
    tProcessed = toc(mainAgentData.tStart);
    printElapsedTime(tProcessed);
    
    % Open enabled plots after processing is finished
    SinglePlotting(mainAgentData.satellites, mainAgentData.storeData, mainAgentData.obs, mainAgentData.settings)    % open enabled plots
    
    % end of epoch-by-epoch processing
    % *************************************************************************
    
    % processing is finished:
    % update plot-panel of GUI_PPP
    handles = disable_plot_checkboxes(handles, mainAgentData.settings);
    path_data4plot = [mainAgentData.settings.PROC.output_dir, '/data4plot.mat'];
    set(handles.edit_plot_path,'String', path_data4plot);
    handles.paths.plotfile = path_data4plot;      % save path to data4plot.mat into handles
    set(handles.pushbutton_load_pos_true,'Enable','On');
    set(handles.pushbutton_load_true_kinematic,'Enable','On');
    handles.paths.lastproc = mainAgentData.settings.PROC.output_dir;     	% save path to last processing into handles
    fprintf('\n---------------------------------------------------------------------\n');
    
    %%%%%%%%%%%%%%%%%% OLD CODE TO BE CLEANED
    old_code =0;
    if old_code
        fprintf('Starting batch processing....\n\n')
        if bool_parfor          % use parfor loop for batch processing
            WaitMessage = parfor_wait(n, 'Waitbar', true);
            
            parfor ii = 1:n                 % loop over rows, process each row
                ROW = TABLE(ii,:);          % current row
                try     % to continue with next file when processing of current fails
                    settings_now = BatchProcessingPreparation(settings, ROW);	% prepare for PPP_main.m
                    % check if mainAgentData.settings for current processing/row are valid
                    valid_settings = checkProcessingSettings(settings_now, false);
                    if ~valid_settings
                        continue            % ... with next processing/row
                    end
                    [~,file,ext] = fileparts(settings_now.INPUT.file_obs);
                    fprintf('\n---------------------------------------------------------------------\n');
                    fprintf('%s%03d%s%03d%s\n%s','Batch Processing File #', ii, ' of ', n, ': ', [file ext]);
                    fprintf('\n---------------------------------------------------------------------\n');
                    % -+-+-+- CALL MAIN FUNCTION  -+-+-+-
                    settings_ = PPP_main(settings_now);         % start processing
                    success(ii) = true;
                    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                end
                WaitMessage.Send;       % update batch processing waitbar
            end
            
            % parfor batch processing is finished
            poolobj = gcp('nocreate');
            delete(poolobj);     % close parallel workers
            WaitMessage.Destroy;        % close batch processing waitbar
            
        else                    % use normal loop for batch processing
            for i = 1:n              % loop over rows, process each row
                ROW = TABLE(i,:);       % current row
                if STOP_CALC
                    break               % -> stop batch processing
                end
                settings_now = BatchProcessingPreparation(settings, ROW);   % prepare for PPP_main.m
                % check if mainAgentData.settings for current processing/row are valid
                valid_settings = checkProcessingSettings(settings_now, false);
                
                %
                if ~valid_settings
                    continue    % with next processing/row
                end
                [~,file,ext] = fileparts(settings_now.INPUT.file_obs);
                fprintf('\n---------------------------------------------------------------------\n');
                fprintf('%s%03d%s%03d%s\n%s\n\n','Batch Processing File #', i, ' of ', n, ': ', [file ext])
                fprintf('\n---------------------------------------------------------------------\n');
                try         % to continue with next file when processing of current fails
                    
                    
                    % -+-+-+- CALL MAIN FUNCTION  -+-+-+-
                    settings_ = PPP_main(settings_now);         % start processing
                    success(i) = true;
                    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                catch ERR_info
                    mess1 = [[file ext] ' failed!'];
                    mess2 = ['Function: ' ERR_info.stack(1).name '.m, line: ' num2str(ERR_info.stack(1).line)];
                    error_str = {mess1, mess2};
                    errordlg(error_str, ['Batch Processing #' sprintf('%02d', i)]);
                end
            end
        end
        
        % Batch processing is over
        l_end = ['End:  ' char(datetime('now'))];
        if ~STOP_CALC   % Batch Processing has finished properly
            failed = TABLE(~success,2);         % check which files have failed
            if isempty(failed)          % everything correctly processed
                mess_header = {'Batch-Processing is done.'; l_start; l_end};
                msgbox(mess_header, 'Achievement', 'help');
                delete('PROCESSLIST/LastFailed.mat');
            else
                mess_header = {'Batch-Processing is done.'; l_start; l_end; 'Failed:'};
                % add number to failed file name (could be vectorized somehow)
                idx_failed = find(~success);
                n_failed = sum(~success);
                failed2 = cell(n_failed, 1);
                for i = 1:n_failed
                    failed2{i} = ['#' num2str(idx_failed(i)) ' ' failed{i}];
                end
                % print messagebox with information about failed files
                msgbox(vertcat(mess_header, failed2), 'Achievement', 'help');
                % save failed processlist in folder PROCESSLIST
                process_list = TABLE(~success, :);
                save('PROCESSLIST/LastFailed.mat', 'process_list');
            end
        else            % Batch Processing was stopped
            [~,file,ext] = fileparts(settings_now.INPUT.file_obs);
            line1 = 'Batch-Processing stopped during:';
            line2 = ['File #' sprintf('%02d', i-1) ', ' file ext];
            info = {line1, line2, l_start, l_end};
            msgbox(info, 'Achievement', 'help')
        end
    end
    
else        % Start Processing of single file
    
    settings = getSettingsFromGUI(handles);         % get input from GUI and save it in structure "settings"
    save([pwd,   '/', 'settings.mat'], 'settings') 	% save mainAgentData.settings from GUI as .mat
    
    % check if mainAgentData.settings for processing are valid
    valid_settings = checkProcessingSettings(settings, false);
    if ~valid_settings
        return
    end
    settings = manipulateProcessingName(settings);
    [~,file,ext] = fileparts(settings.INPUT.file_obs);
    fprintf('\n---------------------------------------------------------------------\n');
    fprintf('%s%s\n\n','Observation file: ',[file ext]);
    
    
    % -+-+-+- CALL MAIN FUNCTION  -+-+-+-
    settings_ = PPP_main(settings);         % start processing
    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    
    % processing is finished:
    % update plot-panel of GUI_PPP
    handles = disable_plot_checkboxes(handles, settings_);
    path_data4plot = [settings_.PROC.output_dir, '/data4plot.mat'];
    set(handles.edit_plot_path,'String', path_data4plot);  	% ||| ugly
    handles.paths.plotfile = path_data4plot;      % save path to data4plot.mat into handles
    set(handles.pushbutton_load_pos_true,'Enable','On');
    set(handles.pushbutton_load_true_kinematic,'Enable','On');
    handles.paths.lastproc = settings_.PROC.output_dir;     	% save path to last processing into handles
    fprintf('\n---------------------------------------------------------------------\n');
end
