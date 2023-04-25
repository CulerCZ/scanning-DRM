[misOriMap_top_kernel_02, oriEBSD_top, misOriAngle_top_02] = errAnalyzeAl( ...
    ebsd_top, indexResult_kernel_1.eulerMap, movingPoints_top, refPoints_top);
[misOriMap_bot_kernel_02, oriEBSD_bot, misOriAngle_bot_02] = errAnalyzeAl( ...
    ebsd_bot, indexResult_kernel_1.eulerMap, movingPoints_bot, refPoints_bot);
%% 
figure, imshow(misOriMap_top_kernel_02,Border="tight")
colormap(jet)
clim([0 20])

figure, imshow(misOriMap_bot_kernel_02,Border="tight")
colormap(jet)
clim([0 20])

misOri = [reshape(misOriMap_top_kernel_02,[],1); reshape(misOriMap_bot_kernel_02,[],1)];
misOri = misOri(~isnan(misOri));
figure,
histogram(misOri,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('prediction error (deg)')
xlim([1 62])
ylim([0 1.5e5])
ylabel('number of pixels')
title('histogram of orientation error')

%% plot error on inverse pole figure
figure, plotIPDF([oriEBSD_top,oriEBSD_bot],[misOriAngle_top_02;misOriAngle_bot_02],vector3d.Z,'points',1e6,'MarkerSize',1)
colormap("jet")

%% functions used in this script
function [misOriMap, oriEBSD, misOriAngle] = errAnalyzeAl( ...
    ebsd, eumap, movingPoints, refPoints, options)
    arguments
        ebsd
        eumap
        movingPoints
        refPoints
        options.GBthres (1,1) double = 10
        options.isPlot (1,1) logical = false
    end
    ebsd_temp = ebsd;
    [grains, ebsd_temp.grainId] = calcGrains(ebsd_temp('indexed'),options.GBthres*degree);
    
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
    % color_ebsd_temp = plot_ipf_map(EUmap_temp);
    
    validRange = [floor(sizeEBSD(1)/4),floor(sizeEBSD(1)/4*3)];
    shiftVal = (find(isIndexedMap(:,validRange(1)),1) - find(isIndexedMap(:,validRange(2)),1)) / (validRange(2)-validRange(1));
    tformMat = [1 0 0; shiftVal 1 0; 0 0 1];
    tform = affinetform2d(tformMat);
    EUmap_ebsd = imwarp(EUmap_temp,tform,'nearest');
    isIndexedMap_ebsd = imwarp(isIndexedMap,tform,"nearest");
    EUmap_ebsd = permute(EUmap_ebsd,[3 1 2]);
    EUmap_ebsd(:,~isIndexedMap_ebsd) = NaN;
    EUmap_ebsd = permute(EUmap_ebsd,[2 3 1]);
    
    % colorDRMoriginal = plot_ipf_map(eumap);
    colorEBSDoriginal = plot_ipf_map(EUmap_ebsd);
    tform_register = fitgeotrans(movingPoints,refPoints,'affine');
    output_region = imref2d(size(colorEBSDoriginal));
    EUmap_trans = imwarp(eumap,tform_register,'nearest','OutputView',output_region);
    if options.isPlot
        figure, imshowpair(plot_ipf_map(EUmap_trans),colorEBSDoriginal,'montage')
    end
    
    [n1,n2,~] = size(EUmap_ebsd);
    eulerDRM = reshape(EUmap_trans,n1*n2,3);
    eulerEBSD = reshape(EUmap_ebsd,n1*n2,3);
    cs = ebsd_temp.CSList{2};
    oriDRM = orientation.byEuler(eulerDRM.*degree,cs);
    oriEBSD = orientation.byEuler(eulerEBSD.*degree,cs);
    rot = rotation.byAxisAngle(vector3d.Y,0*degree);
    oriEBSD = rot*oriEBSD;
    rot = rotation.byAxisAngle(vector3d.Z,90*degree);
    oriEBSD = rot*oriEBSD;
    
    misOriAngle = angle(oriDRM,oriEBSD,cs)./degree;
    misOriMap = reshape(misOriAngle,n1,n2);
    misOriMap(~isIndexedMap_ebsd) = NaN;

    if options.isPlot
        % plot indexing error mapping
        figure, imshow(misOriMap,Border="tight")
        colormap(jet)
        clim([0 20])
        colorbar
        % plot pixel-wise indexing error histogram
        figure, histogram(misOriMap,61,'BinLimits',[1 62],...
            'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
        set(gca,'LineWidth',1.5,'FontSize',14)
        xlabel('prediction error (deg)')
        xlim([1 62])
        ylim([0 7e4])
        ylabel('number of pixels')
        title('histogram of orientation error')
    end

    grainSizeThres = 20;
    num_grain = max(grainIdMap,[],'all');
    misOriGrain = zeros(num_grain,1);
    for ii = 1:num_grain
        if sum(grainIdMap == ii,"all") < grainSizeThres
            misOriGrain(ii) = NaN;
            continue
        else
            misOri_temp = misOriMap(grainIdMap == ii);
            misOriGrain(ii) = median(misOri_temp(misOri_temp < prctile(misOri_temp,80)));
        end
    end
    % 
    % figure, histogram(misOriGrain,61,'BinLimits',[1 62],...
    %     'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
    % set(gca,'LineWidth',1.5,'FontSize',14)
    % xlabel('prediction error (deg)')
    % xlim([1 62])
    % ylabel('number of grains')
    % title('histogram of orientation error')

    % figure, plotIPDF(oriEBSD,misOriAngle,vector3d.Z,'points',length(oriEBSD),'MarkerSize',1)
    % colormap("jet")
end