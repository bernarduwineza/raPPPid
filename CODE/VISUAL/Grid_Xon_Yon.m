function [] = Grid_Xon_Yon()
    % function to turn both the x-grid and the y-grid on for the current plot
    %
    % This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner
    % *************************************************************************

    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
end