

% Fit planes onto epochwise VTEC estimates
n = 1; 
plv_gim = [lat_pp(1:n,1:32)*180/pi; lon_pp(1:n,1:32)*180/pi; vtec_gim(1:n,1:32)];
plv_phi = [lat_pp(1:n,1:32)*180/pi; lon_pp(1:n,1:32)*180/pi; vtec_phi(1:n,1:32)];

B = plv_gim'; 
B = rmmissing(B);
s_fit = fit([B(:,1), B(:,2)], B(:,3),'poly13'); 
a = plot(s_fit,[B(:,1), B(:,2)], B(:,3)); 
h = colorbar;
h.Label.String = "VTEC (TECu)";
h.Label.FontSize = 16;
view([90, -90])
colormap(roma100)
hold on 
% scatter3(B(:,1), B(:,2), B(:,3), 200, 'red', 'filled', 'Marker','o')
xlabel('Lat (deg)');
ylabel('Lon (deg)');
zlabel('VTEC (TECu)');
%%
B = plv_phi'; 
B = rmmissing(B); 
s_fit = fit([B(:,1), B(:,2)], B(:,3),'poly13'); 
a = plot(s_fit,[B(:,1), B(:,2)], B(:,3)); 
h = colorbar;
h.Label.String = "VTEC (TECu)";
h.Label.FontSize = 16;
view([90, -90])
colormap(roma100)
xlabel('Lat (deg)');
ylabel('Lon (deg)');
zlabel('VTEC (TECu)');
%%
[n_gim,~,p_gim] = affine_fit(plv_gim');
[n_phi,~,p_phi] = affine_fit(plv_phi');

fig_title = 'GIM Iono Estimates: Plane fit';
togglefig(fig_title); clf;

plot3(plv_gim(1,:), plv_gim(2,:), plv_gim(3,:),'r.', 'MarkerSize',15, 'MarkerFaceColor','red');
legend('GIM (VTEC)', 'Phase Iono (VTEC)')
hold on;
plot3(plv_phi(1,:), plv_phi(2,:), plv_phi(3,:),'b.', 'MarkerSize',15, 'MarkerFaceColor','blue'); 

%plot the two points p_1 and p_2
plot3(p_gim(1),p_gim(2),p_gim(3),'ro','markersize',5,'markerfacecolor','red');
plot3(p_phi(1),p_phi(2),p_phi(3),'bo','markersize',5,'markerfacecolor','blue');


% %plot the normal vector
quiver3(p_gim(1),p_gim(2),p_gim(3),n_gim(1)/3,n_gim(2)/3,n_gim(3)/3,'r','linewidth',2)
h = quiver3(p_phi(1),p_phi(2),p_phi(3),n_phi(1)/3,n_phi(2)/3,n_phi(3)/3,'b','linewidth',2);

%first plane
[X,Y] = meshgrid(plv_gim(1,:), plv_gim(2,:));
s = surf(X,Y, -(n_gim(1)/n_gim(3)*X + n_gim(2)/n_gim(3)*Y - dot(n_gim,p_gim)/n_gim(3)),'facecolor','red','facealpha',0.5);
s.EdgeColor = 'none';
s = surf(X,Y, - (n_phi(1)/n_phi(3)*X + n_phi(2)/n_phi(3)*Y - dot(n_phi,p_phi)/n_phi(3)),'facecolor','blue','facealpha',0.5);
s.EdgeColor = 'none';

grid minor
zlim([0 15])

legend('GIM (VTEC)', 'Phase Iono (VTEC)')
xlabel('Lat (deg)');
ylabel('Lon (deg)');
zlabel('VTEC (TECu)');

angle = acosd(dot(n_gim,n_phi));
if angle>90
    angle = 180-angle;
end
disp(angle) 


%% 
B = plv_gim'; 
B = rmmissing(B); 
s_fit = fit([B(:,1), B(:,2)], B(:,3),'poly23'); 
plot(s_fit,[B(:,1), B(:,2)], B(:,3)); 
hold on 
