function input = get_V3GR(obs, input, settings)
% Defines and reads all VMF3+GRAD information, collected under the term V3GR.
%
% Coded:
%   01 Feb 2019 by D. Landskron
%
% Revision:
%   ...
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************



% define which jd's are needed
V3GR_data_all = [];
jd_all = (floor(obs.startdate_jd-0.5)+0.5-1 : floor(obs.startdate_jd-0.5)+0.5+1);   % get the 3 jd's: previous day, current day and subsequent day (a bit complicated because a julian day starts at *.5)

% find out, which products are to be read
if ~strcmpi(settings.TROPO.Gh,'GRAD')   &&   ~strcmpi(settings.TROPO.Gw,'GRAD')
    product = 'VMF3';
elseif ~strcmpi(settings.TROPO.zhd,'VMF3')   &&   ~strcmpi(settings.TROPO.zwd,'VMF3')   &&   ~strcmpi(settings.TROPO.mfh,'VMF3')   &&   ~strcmpi(settings.TROPO.mfw,'VMF3')
    product = 'GRAD';
else
    product = 'VMF3+GRAD';
end

% start the loop over all necessary daily files
for i_jd = 1:length(jd_all)
    
    % first determine all date variables
    mjd = jd_all(i_jd) - 2400000.5;
    [ year, month, day] = jd2cal_GT(jd_all(i_jd));
    year_str = num2str(year);
    doy = jd2doy_GT(jd_all(i_jd));
    doy_str = num2str(doy,'%03d');
    
    % define the subdirectory
    if year < 2008
        NWM_source = 'EI';
    else
        NWM_source = 'OP';
    end
    
    % define the paths for V3GR
    dir_V3GR_sitewise = ['../DATA/V3GR/sitewise/V3GR_' NWM_source '/'];
    url_V3GR_sitewise = ['https://vmf.geo.tuwien.ac.at/trop_products/GNSS/V3GR/V3GR_' NWM_source '/daily/'];
    dir_V3GR_gridwise = ['../DATA/V3GR/gridwise/V3GR_' NWM_source '/'];
    url_V3GR_gridwise = ['https://vmf.geo.tuwien.ac.at/trop_products/GRID/5x5/V3GR/V3GR_' NWM_source '/'];
    
    
    file_V3GR = [year_str doy_str '.v3gr_g'];
    
    % check if the respective V3GR yearly file(s) is locally available. If not, then download. Also update an outdated existing yearly file.
    if ~exist([dir_V3GR_sitewise '/' year_str '/' file_V3GR],'file')   % if file does not exists, download in any case
        mkdir([dir_V3GR_sitewise '/' year_str '/']);
        urlwrite([url_V3GR_sitewise '/' year_str '/' file_V3GR], [dir_V3GR_sitewise '/' year_str '/' file_V3GR]);
        if ~settings.INPUT.bool_parfor
            fprintf('  %s%s%s%s\n' , file_V3GR,' successfully downloaded into ',dir_V3GR_sitewise,'.');
        end
    end
    
    
    % open the V3GR file
    fid = fopen([dir_V3GR_sitewise '/' year_str '/' file_V3GR]);
    V3GR_data = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f','CommentStyle','#');
    fclose(fid);
    
    % concatenate the files, if necessary
    V3GR_data_all = [V3GR_data_all;V3GR_data];
    
end


% if more than one year, then concatenate
clear V3GR_data
for i_col = 1:size(V3GR_data_all,2)
    V3GR_data{i_col} = cat(1, V3GR_data_all{:,i_col});
end

% reduce to the IGS station; if no IGS station, then use grid-wise V3GR instead
ind = ismember(V3GR_data{1,1}, strtrim(obs.stationname));
if sum(ind)~=0
    
    % write to command window that the sitewise V3GR is used
    input.TROPO.V3GR.version = 'sitewise';
    if ~settings.INPUT.bool_parfor
        fprintf('  %s%s%s%s%s%s%s\n','Station ',upper(obs.stationname),' is contained in the ',product,' station list => sitewise ',product,' coefficients applied');
    end
    
    for i_col = 1:length(V3GR_data)
        input.TROPO.V3GR.data{i_col} = V3GR_data{i_col}(ind);
    end
    
else
    
    % write to command window that the gridwise V3GR is used
    input.TROPO.V3GR.version = 'gridwise';
    if ~settings.INPUT.bool_parfor
        fprintf('  %s%s%s%s%s%s%s\n','Station ',upper(obs.stationname),' is not contained in the ',product,' station list => gridwise ',product,' coefficients applied');
    end
    
    % determine all necessary mjd's
    input.TROPO.V3GR.data{2} = unique(V3GR_data{2});
    
    % get the correct VMF3 and/or GRAD coefficients
    V3GR_grid_file = []; 
    approx_pos_WGS84 = cart2geo(settings.INPUT.pos_approx);
    for i_jd = 1:length(input.TROPO.V3GR.data{2})
        
        input.TROPO.V3GR.data{1}{i_jd} = upper(obs.stationname);
        [ input.TROPO.V3GR.data{3}(i_jd,1), input.TROPO.V3GR.data{4}(i_jd,1), input.TROPO.V3GR.data{5}(i_jd,1), input.TROPO.V3GR.data{6}(i_jd,1), ...
          input.TROPO.V3GR.data{10}(i_jd,1), input.TROPO.V3GR.data{11}(i_jd,1), input.TROPO.V3GR.data{12}(i_jd,1), input.TROPO.V3GR.data{13}(i_jd,1), V3GR_grid_file ] = ...
          v3gr_grid_adapted ( dir_V3GR_gridwise, url_V3GR_gridwise, V3GR_grid_file, input.TROPO.V3GR.data{2}(i_jd), approx_pos_WGS84.ph, approx_pos_WGS84.la, approx_pos_WGS84.h, 5 );
        
    end
    
end