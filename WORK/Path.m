classdef Path
% Class for paths which are used in VieVS PPP. Any changes of paths
% have to be done here (only). This file should be in the same folder as
% vievs.m which is used for starting (the GUI of) VieVS PPP.
%  
%   Revision:
%   ...
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
% *************************************************************************
    
    
    properties (Constant = true)
        

        % Programm-folder
        CODE = '../CODE/'
        
        % Directory of data
        DATA     = '../DATA/';
        
        % Directory of results
        RESULTS = '../RESULTS/';
        
        % Path to folder of 7-zip.exe
        ZIP7 = '../CODE/7ZIP/7za.exe';
        
        % Path of TUW ionosphere models
        TUW_IONO = 'A:/Datapool/VieVS-IONO/IONEX/';

    end
    
    methods (Access = public)
    end
end