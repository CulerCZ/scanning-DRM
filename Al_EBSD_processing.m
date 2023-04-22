%% EBSD stitching
rot = rotation.byAxisAngle(vector3d.Y,0*degree);
ebsd_temp = rotate(ebsd_top,rot);
[grains, ebsd_temp.grainId] = calcGrains(ebsd_temp('indexed'),10*degree);
figure, plot(ebsd_temp('indexed'),ebsd_temp('indexed').orientations,'micronbar','off');
hold on
plot(grains.boundary,'LineWidth',2);
hold off
% original ebsd plot of the sample
%%
stepSize = max(ebsd_temp.unitCell) - min(ebsd_temp.unitCell);
pos_x = fix(ebsd_temp.prop.x);
pos_y = fix(ebsd_temp.prop.y);
pos_x = pos_x - min(pos_x);
pos_y = pos_y - min(pos_y);
idxEBSD = [pos_x, pos_y] ./ stepSize + [1,1];
% under the current input settings, a 90-deg clockwise rotation is needed.
sizeEBSD = max(idxEBSD);
clear EUmap_temp
EUmap_temp(1,:,:) = rot90(reshape(ebsd_temp.rotations.phi1,sizeEBSD(1),sizeEBSD(2)),1)./degree;
EUmap_temp(2,:,:) = rot90(reshape(ebsd_temp.rotations.Phi,sizeEBSD(1),sizeEBSD(2)),1)./degree;
EUmap_temp(3,:,:) = rot90(reshape(ebsd_temp.rotations.phi2,sizeEBSD(1),sizeEBSD(2)),1)./degree;

isIndexedMap = rot90(reshape(ebsd_temp.isIndexed,sizeEBSD(1),sizeEBSD(2)),1);
grainIdMap = rot90(reshape(ebsd_temp.grainId,sizeEBSD(1),sizeEBSD(2)),1);

EUmap_temp(:,~isIndexedMap) = NaN;
EUmap_temp = permute(EUmap_temp,[2 3 1]);
color_ebsd_temp = plot_ipf_map(EUmap_temp);
figure, imshow(color_ebsd_temp,'border','tight')
%% ----------------------------------------------------------------------------
% this part of transformation needs careful consideration
validRange = [floor(sizeEBSD(1)/4),floor(sizeEBSD(1)/4*3)];
shiftVal = (find(isIndexedMap(:,validRange(1)),1) - find(isIndexedMap(:,validRange(2)),1)) / (validRange(2)-validRange(1));
tformMat = [1 0 0; shiftVal 1 0; 0 0 1];
tform = affinetform2d(tformMat);
EUmap_ebsd = imwarp(EUmap_temp,tform,'nearest');
isIndexedMap_ebsd = imwarp(isIndexedMap,tform,"nearest");
EUmap_ebsd = permute(EUmap_ebsd,[3 1 2]);
EUmap_ebsd(:,~isIndexedMap_ebsd) = NaN;
EUmap_ebsd = permute(EUmap_ebsd,[2 3 1]);
figure, imshow(imwarp(color_ebsd_temp,tform,"nearest"),'Border','tight');
% ----------------------------------------------------------------------------
% 
% Scale bar to be added


%% register EBSD and DRM dataset and compare indexing error
% register two datasets
EUmap = eumap;
colorDRMoriginal = plot_ipf_map(eumap);
colorEBSDoriginal = plot_ipf_map(EUmap_ebsd);
if ~exist("movingPoints",'var')
    [movingPoints, refPoints] = cpselect(colorDRMoriginal,colorEBSDoriginal,'Wait',true);
end
tform_register = fitgeotrans(movingPoints_top,refPoints_top,'affine');
output_region = imref2d(size(colorEBSDoriginal));
EUmap_trans = imwarp(eumap,tform_register,'nearest','OutputView',output_region);
figure, imshowpair(plot_ipf_map(EUmap_trans),colorEBSDoriginal,'montage')

%% compare EUmap_trans and EUmap_ebsd
[n1,n2,~] = size(EUmap_ebsd);
eulerDRM = reshape(EUmap_trans,n1*n2,3);
eulerEBSD = reshape(EUmap_ebsd,n1*n2,3);
cs = ebsd_temp.CSList{2};
oriDRM = orientation.byEuler(eulerDRM.*degree,cs);
oriEBSD = orientation.byEuler(eulerEBSD.*degree,cs);
rot = rotation.byAxisAngle(vector3d.Y,180*degree);
oriEBSD = rot*oriEBSD;
rot = rotation.byAxisAngle(vector3d.Z,180*degree);
oriEBSD = rot*oriEBSD;

misOriAngle = angle(oriDRM,oriEBSD,cs)./degree;
misOriMap = reshape(misOriAngle,n1,n2);
misOriMap(~isIndexedMap_ebsd) = NaN;
% plot indexing error mapping
figure, imshow(misOriMap,Border="tight")
colormap(jet)
clim([0 20])
colorbar
% plot pixel-wise indexing error histogram
figure,
histogram(misOriMap,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('prediction error (deg)')
xlim([1 62])
ylabel('number of pixels')
title('histogram of orientation error')
%% calculate grain-wise error and plot corresponding histogram
grainSizeThres = 20;
num_grain = max(grainIdMap,[],'all');
misOriGrain = zeros(num_grain,1);
for ii = 1:num_grain
    if sum(grainIdMap == ii,"all") < grainSizeThres
        misOriGrain(ii) = NaN;
        continue
    else
        misOri_temp = misOriMap(grainIdMap == ii);
        misOriGrain(ii) = mean(misOri_temp(misOri_temp < prctile(misOri_temp,90)));
    end
end
figure,
histogram(misOriGrain,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('prediction error (deg)')
xlim([1 62])
ylabel('number of pixels')
title('histogram of orientation error')