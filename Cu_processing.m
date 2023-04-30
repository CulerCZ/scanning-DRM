%% scanning DRM on EJM-Cu sample
% Load an image
figure,
fig_temp = mean(dataNorm.drpMap,3);
imshow(fig_temp,[min(fig_temp,[],"all"),max(fig_temp,[],"all")],'Border','tight');

% Use ginput to get the coordinates of a point
[x, y] = ginput(2);

% Display the coordinates in the MATLAB command window
fprintf('Pixel coordinates: (%.1f, %.1f)\n', x(1), y(1));
fprintf('Pixel coordinates: (%.1f, %.1f)\n', x(2), y(2));

%% line scan over the sample
xlim = [30 390];
ylim = [35 215];
x_sample = fix(linspace(xlim(1),xlim(2),10));
y_sample = fix(linspace(ylim(1),ylim(2),10));
% plot sampling on the raw sample image
figure,
imshow(fig_temp,[min(fig_temp,[],"all"),max(fig_temp,[],"all")],'Border','tight');
[x_temp, y_temp] = meshgrid(y_sample(2:2:10), x_sample(2:2:10));
hold on
scatter(x_temp, y_temp, 'r', 'filled')
hold off

% plot sample DRPs
figure(Position=[0 0 1000 1000])
tiledlayout(5,5,"TileSpacing","compact",Padding="compact")
for ii = 1:5
    for jj = 1:5
        nexttile(5*(ii-1)+jj)
        plotDRP(dataNorm.drpMap(x_sample(2*ii),y_sample(2*jj),:),posInfo);
    end
end


%% line profile over the EJM scanning direction
dataPoints = zeros(length(xlim(1):xlim(2)),5);
for x_temp = xlim(1):xlim(2)
    drp_data = dataNorm.drpMap(x_temp,ylim(1):ylim(2),:);
    for jj = 75:5:95
        dataPoints(x_temp-xlim(1)+1,(jj-75)/5+1) = prctile(drp_data,jj,"all");
    end
    dataPoints(x_temp-xlim(1)+1,6) = mean(drp_data,"all");
end
color = {"#8dd3c7"; "#6dbcc3"; "#4d9eb6"; "#3d81a2"; "#2c6282"};
figure(Position=[0 0 1000 500]),
hold on
for ii = 1:5
    plot(((xlim(1):xlim(2))-xlim(1))*50/1000,dataPoints(:,ii),"LineWidth",2,"Color",color{ii})
end
% plot(((xlim(1):xlim(2))-xlim(1))*50/1000,dataPoints(:,6),"LineWidth",2,"Color","black")
box on
legend(["75%","80%","85%","90%","95%"])
set(gca,"LineWidth",2,"FontSize",14,'ylim',[-0.05 1.05])
%% #8dd3c7, #6dbcc3, #4d9eb6, #3d81a2, #2c6282
