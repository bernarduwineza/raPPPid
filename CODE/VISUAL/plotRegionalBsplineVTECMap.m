% Plot the VTEC map using estimated Bspline coefficients
% (This assumes that a solution from a previous run has bee loaded)

% Define the lat, lon range for the estimated values, i.e. define local
% region.
step_sz = 1;

lat_min = settings.IONO.Bspline.lat_min+1; 
lat_max = settings.IONO.Bspline.lat_max;
lon_min = settings.IONO.Bspline.lon_min;
lon_max = settings.IONO.Bspline.lon_max;

lat_range = lat_min:step_sz:lat_max;
lon_range = lon_min:step_sz:lon_max;

% number of basis functions
kj1 = settings.IONO.Bspline.kj1;
kj2 = settings.IONO.Bspline.kj2;

% Resolution levels
j1 = settings.IONO.Bspline.lat_res_level;
j2 = settings.IONO.Bspline.lon_res_level;

% Degree of the B-spline basis functions
lat_degree = settings.IONO.Bspline.lat_degree;
lon_degree = settings.IONO.Bspline.lon_degree;

% Form knot points (using Bspline toolbox )
lat_knots = augknt(linspace(lat_min, lat_max, kj1-lat_degree+1), j1);
lon_knots = augknt(linspace(lon_min, lon_max, kj2-lon_degree+1), j2);

% Get global VTEC from IONEX file
path_ionex = settings.IONO.file_ionex;
ionex = read_ionex_TUW(settings.IONO.file_ionex);

colormap jet;
fprintf('Epoch 000000');
clf;
for i =1:length(storeData.bspline_coeff)
    if mod(i,15)==0
        idx_sats = storeData.lat_pp(i,:) ~= 0;

        lat_pp = storeData.lat_pp(i,idx_sats); 
        lon_pp = storeData.lon_pp(i,idx_sats); 

        lat_ctrl_pts = lat_min:5:lat_max; 
        lon_ctrl_pts = lon_min:5:lon_max;

        Ttr = model_save.Ttr(i,3);       	% time of emission [sow (seconds of week)] = time of obs. - runtime
        
        % interpolate VTEC
        count = 0; 
        [B,A] = ndgrid(lon_ctrl_pts, lat_ctrl_pts);
        ctrl_pts = [A(:),B(:)];
        vtec_ctrl_pt_gim = zeros(length(ctrl_pts),1); 
        for ii = 1:length(ctrl_pts)
            vtec_ctrl_pt_gim(ii) = iono_gims(ctrl_pts(ii,1), ctrl_pts(ii,2), Ttr, ionex, settings.IONO.interpol);
        end

        lat_lon_bases = compute_bspline_bases(settings, lat_pp*pi/180, lon_pp*pi/180);
        [lat_bases, lon_bases] = compute_lat_lon_bspline_bases(settings, lat_range'*pi/180, lon_range'*pi/180);
        
        lat_lon_bases_ctrl_pts = compute_bspline_bases(settings, ctrl_pts(:,1)*pi/180, ctrl_pts(:,2)*pi/180);
        lat_lon_bases_ctrl_pts(isnan(lat_lon_bases_ctrl_pts)) = 0; 

        % Estimate Bspline coefficients using LS
        vtec = storeData.iono_vtec(i, storeData.iono_vtec(i,:)~=0)'; 
        % coef_iono_LS = lat_lon_bases'\vtec;
        % coef_iono_LS = reshape(coef_iono_LS, [kj1, kj2]);
        % coef_iono_ls = reshape(lsqr(lat_lon_bases_ctrl_pts',vtec_ctrl_pt_gim, 1e-6), [kj1, kj2]);
        coef_iono_ls    = reshape(lat_lon_bases_ctrl_pts'\vtec_ctrl_pt_gim, [kj1, kj2]);
        % The estimated Bspline coefficients using the general KF 
        coef_iono_kf= storeData.bspline_coeff(i,:);
        coef_iono_kf = reshape(coef_iono_kf', [kj1, kj2]);

        % Compute the VTEC values for every point in the range
        psi_lat = spcol(lat_knots, j1, lat_range); 
        psi_lon = spcol(lon_knots, j2, lon_range);
        % vtec_vals_kf = psi_lat * coef_iono * psi_lon';
        vtec_vals_kf = lat_bases'*coef_iono_kf*lon_bases; 
        vtec_vals_ls = lat_bases'*coef_iono_ls*lon_bases;

        PLOT_TEC = 1; 
        if PLOT_TEC
            % Plot VTEC values
            plt = surf(lat_range,lon_range, vtec_vals_kf.');

            hold on;
            plot3(lat_pp, lon_pp, vtec, 'r.', 'MarkerSize',12);
            hold off;
            clim([0 20]);
            % zlim([0 100])

            xlabel('Lat (deg)', 'FontSize',14)
            ylabel('Lon (deg)', 'FontSize',14)
            cbar = colorbar;
            cbar.Label.String = 'VTEC (TECu)';
            cbar.FontSize = 14;
            view([90 270]);
            % view([-180 10]);
            drawnow limitrate;
            fprintf('\b\b\b\b\b\b%5d\n', i);
        else 
            plot(i,coef_iono_kf(:), 'b.');
            hold on;
        end 
        drawnow;
    end
end
fprintf('\n');
drawnow;

% plt.EdgeColor = 'interp';
% plt.FaceColor = 'interp';