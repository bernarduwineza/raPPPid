%% Plot the KF coordinate estimate and covariances
% 
% It is assumed that the `storeData` variable is available in the workspace.
%
% (c) 2023, J.-B. Uwineza


figure(1); clf; 


% x coordinate 
x_coord = storeData.param(:,1)';
k = length(x_coord);
x_coord(isnan(x_coord)) = 0;

x_coord = x_coord - mean(x_coord);
x_coord_var = storeData.param_var(:,1)'; 

patches = [x_coord+x_coord_var; x_coord-x_coord_var]';
plot(x_coord, 'b.'); 
hold on; 
plotshaded(linspace(1,k,k), patches, 'r', 0.7);

% ylim([-1 1]);
grid minor
