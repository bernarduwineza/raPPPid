function [lat_bases, lon_bases] = compute_lat_lon_bspline_bases(settings, lat_pp, lon_pp)

kj1 = settings.IONO.Bspline.kj1;                             % number of basis functions
kj2 = settings.IONO.Bspline.kj2;

lat_degree = settings.IONO.Bspline.lat_degree; % Degree of the B-spline basis functions
lon_degree = settings.IONO.Bspline.lon_degree;

lat_min = settings.IONO.Bspline.lat_min; 
lat_max = settings.IONO.Bspline.lat_max;
lon_min = settings.IONO.Bspline.lon_min;
lon_max = settings.IONO.Bspline.lon_max;

lat_knots = linspace(lat_min, lat_max, kj1-lat_degree+1);
lat_knots = [repmat(lat_min, 1, lat_degree), lat_knots, repmat(lat_max, 1, lat_degree)];

lon_knots = linspace(lon_min, lon_max, kj2-lon_degree+1);
lon_knots = [repmat(lon_min, 1, lon_degree), lon_knots, repmat(lon_max, 1, lon_degree)];

% Initialize a matrix to store the basis function values for each basis
lat_bases = zeros(kj1, length(lat_pp));
lon_bases = zeros(kj2, length(lon_pp));

% Calculate B-spline basis functions
for i = 1:kj1
    lat_bases(i, :) = bspline_basis(i, lat_degree, lat_knots, lat_pp(:,1)*180/pi);
end
for i = 1:kj2
    lon_bases(i, :) = bspline_basis(i, lon_degree, lon_knots, lon_pp(:,1)*180/pi);
end
end