%% Plot VTEC map from an ionex file 

ionex = read_ionex_TUW('../../DATA/IONO/2019/280/igsg2800.19i');
load('/Users/Bernard/Downloads/ScientificColourMaps8/roma/DiscretePalettes/roma100.mat')
set(0,'DefaultFigureWindowStyle','docked')
no_maps = size(ionex.map,3);
koeff = 10^ionex.exponent;

lats = ionex.lat(2): -ionex.lat(3): ionex.lat(1);
lons = ionex.lon(1): ionex.lon(3): ionex.lon(2);

lat_int = 5;
lon_int = 4;

for i = 1:no_maps
    figure(i)
    iono_map = ionex.map(:,:,i)*koeff;
    im = imagesc(iono_map);
    % im.Interpolation = 'bilinear';
    clim([0 40])
    % colormap jet
    colormap(flip(roma100))
    colorbar
    xticks(1:lon_int:length(lons));
    xticklabels(lons(1:lon_int:end));

    yticks(1:lat_int:length(lats));
    yticklabels(-lats(1:lat_int:end));

    h = colorbar;
    h.Label.String = "VTEC (TECu)";
    h.Label.FontSize = 16;
    xlabel('Longitude (deg)')
    ylabel('Latitute (deg)')
end