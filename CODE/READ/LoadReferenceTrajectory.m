function [pos_ref_geo, North_ref, East_ref] = ...
    LoadReferenceTrajectory(filepath, leap_sec, gpstime)
% Loads the data from a reference trajectory. 
% 
% INPUT:
%   filepath        string, full relative filepath to reference trajectory
%   leap_sec        integer, number of leap seconds of processing
%   gpstime         vector, gps time (sow) of processing's epochs
% OUTPUT: reference positions interpolated from the reference trajectory to
%         the points in time of the PPP solution
%	pos_ref_geo     struct [.ph, .la, h], WGS84 reference points
%   North_ref       vector, UTM North reference points
%   East_ref        vector, UTM East reference points
%
% Revision:
%   ...
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************


%% Load reference trajectory
[~, ~, ext] = fileparts(filepath);
switch ext
    
    case '.nmea'
        [utc, lat_wgs84, lon_wgs84, h_wgs84] = ReadNMEAFile(filepath);
        sod_true = utc + leap_sec;          % convert utc to gps time
        lat_wgs84 = lat_wgs84 / 180 * pi;   % convert [°] to [rad]
        lon_wgs84 = lon_wgs84 / 180 * pi;   % convert [°] to [rad]
        [North_true, East_true] = ell2utm_GT(lat_wgs84, lon_wgs84);
        
    case '.txt'
        % simple text file with columns: gpstime | X | Y | Z
        fid = fopen(filepath);
        DATA = textscan(fid,'%f %f %f %f');
        fclose(fid);
        % save data
        sow = DATA{:,1}; X = DATA{:,2}; Y = DATA{:,3}; Z = DATA{:,4};
        sod_true = mod(sow, 86400);     % convert sow into sod
        % convert coordinates
        n = numel(sod_true);
        lat_wgs84 = NaN(n,1); lon_wgs84 = lat_wgs84; h_wgs84 = lat_wgs84;
        North_true = lat_wgs84; East_true = lat_wgs84;
        for i = 1:n
            temp_geo = cart2geo([X(i), Y(i), Z(i)]);
            [North_true(i), East_true(i)] = ell2utm_GT(temp_geo.ph, temp_geo.la);
            lat_wgs84(i) = temp_geo.ph;
            lon_wgs84(i) = temp_geo.la;
            h_wgs84(i) = temp_geo.h;
        end
        
    case '.pos'
        % text file with columns:
        % yyyy/mm/dd | hh:mm:ss.sss | lat [°] | lon [°] | height [m] | ignore
        % time is already GPS time
        fid = fopen(filepath);
        header = fgetl(fid);
        
        switch header
            case '%  GPST                   lat(deg)      lon(deg)         h(m)      OK' 
                DATA = textscan(fid,'%f/%f/%f %f:%f:%f %f %f %f %f', 'HeaderLines', 0);
                fclose(fid);
                lat_wgs84 = DATA{:,7} / 180 * pi;   % convert [°] to [rad]
                lon_wgs84 = DATA{:,8} / 180 * pi;   % convert [°] to [rad]
                h_wgs84   = DATA{:,9};
                sod_true = DATA{:,4} * 3600 + DATA{:,5}*60 + DATA{:,6};     % already GPS time
                [North_true, East_true] = ell2utm_GT(lat_wgs84, lon_wgs84);
                
            otherwise
                pos_ref_geo.ph = NaN;
                pos_ref_geo.la = NaN;
                pos_ref_geo.h  = NaN;
                North_ref = NaN;
                East_ref = NaN;
                errordlg('Implement read-in function! Check LoadReferenceTrajectory.', 'Error');
                return
        end
        
    otherwise
        pos_ref_geo.ph = NaN;
        pos_ref_geo.la = NaN;
        pos_ref_geo.h  = NaN;
        North_ref = NaN;
        East_ref = NaN;
        errordlg('Implement read-in function! Check LoadReferenceTrajectory.', 'Error');
        return
        
end



%% interpolate reference trajectory to the time-stamps of the PPP solution
sod = mod(gpstime, 86400);
% creates timeseries variables for resampling
tseries_lat =   timeseries(lat_wgs84,  sod_true);
tseries_lon =   timeseries(lon_wgs84,  sod_true);
tseries_h   =   timeseries(h_wgs84,    sod_true);
tseries_North = timeseries(North_true, sod_true);
tseries_East =  timeseries(East_true,  sod_true);
% resample at the time of the PPP epochs
lat_ref = resample(tseries_lat, sod);
lon_ref = resample(tseries_lon, sod);
h_ref   = resample(tseries_h, sod);
N_ref   = resample(tseries_North, sod);
E_ref   = resample(tseries_East, sod);
% save in output variables
pos_ref_geo.ph = lat_ref.Data;
pos_ref_geo.la = lon_ref.Data;
pos_ref_geo.h  = h_ref.Data;
North_ref      = N_ref.Data;
East_ref       = E_ref.Data;
