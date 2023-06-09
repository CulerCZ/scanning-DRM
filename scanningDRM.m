%% scanning DRM data processing script

[dataRaw, dataBG, posInfo] = loadRawData(...
    AlfromUKRaw20230418181452, AlfromUKBG20230418181452, pixel_coord);
offset = 40;

dataNorm = generateData(dataRaw, dataBG, posInfo, offset=offset);
%%
posInfo_new.phi = pixel_coord(1,:);
posInfo_new.theta = pixel_coord(2,:);
posInfo_new.phi(posInfo_new.theta<=70) = rem(270-posInfo_new.phi(posInfo_new.theta<=70)+360, 360);  % under current settings
posInfo_new.phi(posInfo_new.theta>70) = rem(-90+posInfo_new.phi(posInfo_new.theta>70)+360,360);
plotPara = calcPlotPixels(posInfo,resNum=20);
dataNorm.plotPara = plotPara;
%%
figure, imagesc(unique(dataNorm.x), unique(dataNorm.y),...
    reshape(median(dataNorm.drplist,2),dataNorm.num_x,dataNorm.num_y))

%% create DRP Library and run dictionary indexing
ang_res = 3;   
drpLib = createDRPLib(posInfo, ang_res*degree, faceting=[1,0,0], ...
    fitting_para=[1,0.7,25,4,0.8,8]);
%%
indexResult = IndexEngine_sDRM(dataNorm, drpLib);
figure, imshow(indexResult.distanceMap,[min(indexResult.distance),max(indexResult.distance)])
colormap(jet)
figure, imshow(plot_ipf_map(indexResult.eulerMap),Border="tight")
%%
figure, plotDRP(dataBG.drp,plotPara,caxis=[min(dataBG.drp),max(dataBG.drp)]);
%% plot indexing results
folderpath = ".\results_temp";  % save result image temporarily
offset_values = 0:3:72;
errmisOri_stack = zeros(length(offset_values),n1,n2);
for ii = 1:length(offset_values)
    dataNorm = generateData(dataRaw, dataBG, posInfo, offset=offset_values(ii));
    indexResult = IndexEngine_sDRM(dataNorm, drpLib);
    eumap = indexResult.eulerMap;
    
    colorDRMoriginal = plot_ipf_map(eumap);
    colorEBSDoriginal = plot_ipf_map(EUmap_ebsd);
    if ~exist("movingPoints",'var')
        [movingPoints, refPoints] = cpselect(colorDRMoriginal,colorEBSDoriginal,'Wait',true);
    end
    tform_register = fitgeotrans(movingPoints,refPoints,'affine');
    output_region = imref2d(size(colorEBSDoriginal));
    EUmap_trans = imwarp(eumap,tform_register,'nearest','OutputView',output_region);
    % figure, imshowpair(plot_ipf_map(EUmap_trans),colorEBSDoriginal,'montage')
    
    % compare EUmap_trans and EUmap_ebsd
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
    % plot indexing error mapping
    f1 = figure("name",sprintf("offset value %d ",offset_values(ii)));
    imshow(misOriMap,Border="tight")
    colormap(jet)
    clim([0 20])
    saveas(f1, fullfile(folderpath,sprintf("err_map_offset_%02d.tif",offset_values(ii))));
    close(f1)
    errmisOri_stack(ii,:,:) = misOriMap;
    fprintf("Save error map %02d / %02d offset values.\n",[ii length(offset_values)]);
    % plot pixel-wise indexing error histogram
    f2=figure("name",sprintf("offset %d error histogram",offset_values(ii)));
    histogram(misOriMap,61,'BinLimits',[1 62],...
        'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
    set(gca,'LineWidth',1.5,'FontSize',14)
    xlabel('prediction error (deg)')
    xlim([1 62])
    ylim([0 7e4])
    ylabel('number of pixels')
    saveas(f2,fullfile(folderpath,sprintf("err_histo_offset_%02d.tif",offset_values(ii))));
    close(f2)
    fprintf("Finishing %02d / %02d offset values.\n",[ii length(offset_values)]);
end
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

%%
% eumap = repositionData(indexResult.Euler, dataNorm.x, dataNorm.y);
% figure, imshow(plot_ipf_map(eumap))
figure(Position=[100 100 800 800])
tiledlayout(4,4,"TileSpacing","compact","padding","compact")
rand_idx = randi(size(drpLib.drpList,1),16);
for idx = 1:length(rand_idx)
    nexttile(idx)
    plotDRP(drpLib.drpList(rand_idx(idx),:), plotPara,format="3d")
end
%%
figure(Position=[100 100 800 800])
tiledlayout(4,4,"TileSpacing","compact","padding","compact")
rand_idx = randi(size(dataNorm.drplist,1),16);
for idx = 1:length(rand_idx)
    nexttile(idx)
    plotDRP(dataNorm.drplist(rand_idx(idx),:), plotPara,format="3d")
    colorbar
end

%% select DRPs to be shown
drp_selected = showSampleDRP(dataNorm, plotPara,format="3d");

%% check indexing results
% checkDetail = checkIndexResult(dataNorm, indexResult);
rot = rotation.byAxisAngle(vector3d.Z,90*degree);
ebsd_temp = rotate(ebsd_top,rot);
% idList = [141151 124938 262255 282098 446028];
idList = [290283 655599 493316 450837 173017];
plotDRPfromEBSD(dataNorm.posInfo,ebsd_temp,idList,plotPara);
%% test new indexing engine with 5 cancidates
k_temp = 10;
dataSelected.drplist = checkDetail.drpSelected;
dataSelected.posInfo = dataNorm.posInfo;
dataSelected.num_x = 1;
dataSelected.num_y = 5;
dataSelected.drpMap = reshape(dataSelected.drplist,1,5,608);
indexResultShort = IndexEngine_sDRM(dataSelected, drpLib, K=k_temp);

% figure(Position=[100 100 200*k_temp 1000])
% tiledlayout(5,k_temp,"TileSpacing","compact","Padding","compact")
% for ii = 1:5
%     for jj = 1:k_temp
%         nexttile((ii-1)*k_temp+jj)
%         plotDRP(drpLib.drpList(indexResultShort.Idx(ii,jj),:),dataNorm.posInfo);
%         fprintf("plotting %d / %d subfigures...\n",[(ii-1)*k_temp+jj, 5*k_temp])
%     end
% end

drpEBSDtruth = zeros(5, dataNorm.num_pixel);
distanceEBSD = zeros(5, 2);
idList = [290283 655599 493316 450837 173017];
figure(Position=[100 100 600 1000])
tiledlayout(5,3,"TileSpacing","compact","Padding","compact")
for ii = 1:5
    rotTemp = ebsd_temp(idList(ii)).rotations;
    eulerTemp = [rotTemp.phi1, rotTemp.Phi, rotTemp.phi2] ./ degree;
    drpEBSDtruth(ii,:) = drpSimCone(posInfo, eulerTemp);
    distanceEBSD(ii,1) = norm(drpEBSDtruth(ii,:)-checkDetail.drpSelected(ii,:));
    distanceEBSD(ii,2) = indexResultShort.distance(ii,1);
    nexttile(ii*3-1)
    plotDRP(drpEBSDtruth(ii,:),plotPara,format="2d")
    nexttile(ii*3-2)
    plotDRP(checkDetail.drpSelected(ii,:),plotPara,format="2d")
    nexttile(ii*3)
    plotDRP(drpLib.drpList(indexResultShort.Idx(ii,1),:),plotPara,format="2d")
end



%% quick test of the functions
dataKernel = kernelSmooth(dataNorm,kernelSize=1);
% figure(Position=[100 100 800 400])
% tiledlayout(2,4,"TileSpacing","compact","padding","compact")
% rand_idx = randi(size(dataKernel.drplist,1),4);
% for idx = 1:length(rand_idx)
%     nexttile(idx)
%     plotDRP(dataNorm.drplist(rand_idx(idx),:), posInfo)
%     nexttile(idx+length(rand_idx))
%     plotDRP(dataKernel.drplist(rand_idx(idx),:), posInfo)
% end

indexResult_kernel_1 = IndexEngine_sDRM(dataKernel, drpLib);
% figure, imshow(plot_ipf_map(indexResult_kernel_2.eulerMap),Border="tight")
%%
figure("Position",[100 100 300 300])
plotDRP(drpSimCone(dataNorm.posInfo, squeeze(indexResult.eulerMap(100,100,:))),plotPara)

%% plot BG subtracted DRPs and raw DRPs
figure(Position=[100 100 1600 400])
tiledlayout(2,8,"TileSpacing","compact","padding","compact")
rand_idx = randi(size(dataRaw.drplist,1),8,1);
for idx = 1:length(rand_idx)
    nexttile(idx)
    drptemp = dataRaw.drplist(rand_idx(idx),:) - dataBG.drp;
    plotDRP(drptemp, plotPara, caxis=[prctile(drptemp,20),prctile(drptemp,95)])
    nexttile(idx+length(rand_idx))
    plotDRP(dataRaw.drplist(rand_idx(idx),:), plotPara, caxis= ...
        [prctile(dataRaw.drplist(rand_idx(idx),:),20) prctile(dataRaw.drplist(rand_idx(idx),:),95)])
end
%% function supporting package
% ----------------------------------------------------------------------------
% function to load data
function [dataRaw, dataBG, posInfo] = loadRawData(...
    datasetRaw, datasetBG, pixel_coord)
% OUTPUTS are:
% dataRaw with (1) position x; (2) position y; (3) corresponding DRP lists
% dataBG with (1) offset; (2) gain; (3) background DRP list
% posInfo with (1) azimuth angle phi; (2) polar angle theta
    dataRaw.x = datasetRaw(:,1);
    dataRaw.y = datasetRaw(:,2);
    dataRaw.drplist = datasetRaw(:,3:end);
    dataBG.offset = datasetBG(1);
    dataBG.gain = datasetBG(2);
    dataBG.drp = datasetBG(3:end);
    posInfo.theta = pixel_coord(2,:);
    posInfo.phi = pixel_coord(1,:);
    posInfo.phi(posInfo.theta<=70) = rem(270-posInfo.phi(posInfo.theta<=70)+360, 360);  % under current settings
    posInfo.phi(posInfo.theta>70) = rem(-90+posInfo.phi(posInfo.theta>70)+360,360);
    fprintf("Raw dataset is loaded!\n");
end


% create dataNorm for further processing
function dataNorm = generateData(dataRaw, dataBG, posInfo, options)
    arguments
        dataRaw
        dataBG
        posInfo
        options.offset (1,1) double = 40
        options.gain (1,1) double = 20
    end
    
    pos_x = sort(unique(dataRaw.x),"ascend");
    pos_y = sort(unique(dataRaw.y),"ascend");
    [xx, yy] = meshgrid(pos_x, pos_y);
    num_x = numel(pos_x);
    num_y = numel(pos_y);
    num_pixel = size(dataRaw.drplist,2);
    
    dataNorm.drpMap = nan(num_x, num_y, num_pixel);
    for idx = 1:length(dataRaw.x)
        x_temp = dataRaw.x(idx);
        y_temp = dataRaw.y(idx);
        dataNorm.drpMap(pos_x == x_temp, pos_y == y_temp, :) = ...
            dataRaw.drplist(idx, :) - dataBG.drp;
    end
    
    dataNorm.x = reshape(xx, [], 1);
    dataNorm.y = reshape(yy, [], 1);
    dataNorm.num_x = num_x;
    dataNorm.num_y = num_y;
    dataNorm.num_pixel = num_pixel;
    
    dataNorm.drplist = reshape(dataNorm.drpMap, [], num_pixel);
    dataNorm.drplist = (dataNorm.drplist + options.offset) * options.gain;
    dataNorm.drplist(dataNorm.drplist < 0) = 0;
    % dataNorm.drpMap = (dataNorm.drpMap + options.offset) * options.gain;
    % dataNorm.drpMap(dataNorm.drpMap < 0) = 0;
    dataNorm.drplist = dataNorm.drplist / prctile(dataNorm.drplist,95,"all");
    % dataNorm.drpMap = dataNorm.drpMap / prctile(dataNorm.drplist,95,"all");
    dataNorm.drplist(dataNorm.drplist > 1) = 1;
    % dataNorm.drpMap(dataNorm.drpMap > 1) = 1;
    dataNorm.drpMap = reshape(dataNorm.drplist,num_x,num_y,num_pixel);
    dataNorm.posInfo = posInfo;

    fprintf("Dataset is reday for further processing!\n");
end


% quick plot DRP in polar coordinates
function plotPara = calcPlotPixels(posInfo,options)
    arguments
        posInfo
        options.resNum (1,1) double = 20
    end
    nn = length(posInfo.phi);
    theta_list_original = sort(unique(posInfo.theta),'ascend');
    theta_list = zeros(1,length(theta_list_original)+1);
    for ii = 1:length(theta_list_original)
        if theta_list_original(ii) <= 70
            theta_list(ii) = theta_list_original(ii);
        else
            theta_list(ii) = (theta_list_original(ii) + theta_list_original(ii-1)) / 2;
        end
    end
    theta_list(end) = 90;
    resNum = options.resNum;
    patchPoints = zeros(2*resNum,3,nn);
    for ii = 1:nn
        % datapoints in the lower range
        if posInfo.theta(ii) <= 70
            phi_temp = [posInfo.phi(ii)-5; posInfo.phi(ii)+5; posInfo.phi(ii)+5; posInfo.phi(ii)-5];
            kk = find(theta_list_original == posInfo.theta(ii));
            theta_temp = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
    
            phi = [linspace(phi_temp(1),phi_temp(2),resNum)'; linspace(phi_temp(3),phi_temp(4),resNum)'];
            theta = [linspace(theta_temp(1),theta_temp(2),resNum)'; linspace(theta_temp(3),theta_temp(4),resNum)'];
            [x,y,z] = sph2cart(phi.*degree, theta.*degree, ones(size(phi)));
            patchPoints(:,:,ii) = [x,y,z];
        else
            switch find(theta_list_original == posInfo.theta(ii))
                case 17
                    jj = floor(posInfo.phi(ii) / 22.5);
                    kk = find(theta_list_original == posInfo.theta(ii));
                    phi = [linspace(22.5*jj,22.5*(jj+1),resNum)'; linspace(22.5*(jj+1),22.5*jj,resNum)'];
                    theta = [linspace(theta_list(kk),theta_list(kk),resNum)'; linspace(theta_list(kk+1),theta_list(kk+1),resNum)'];
                    [x,y,z] = sph2cart(phi.*degree, theta.*degree, ones(size(phi)));
                    patchPoints(:,:,ii) = [x,y,z];
                case 18
                    jj = floor(posInfo.phi(ii) / 30);
                    kk = find(theta_list_original == posInfo.theta(ii));
                    phi = [linspace(30*jj,30*(jj+1),resNum)'; linspace(30*(jj+1),30*jj,resNum)'];
                    theta = [linspace(theta_list(kk),theta_list(kk),resNum)'; linspace(theta_list(kk+1),theta_list(kk+1),resNum)'];
                    [x,y,z] = sph2cart(phi.*degree, theta.*degree, ones(size(phi)));
                    patchPoints(:,:,ii) = [x,y,z];
                case 19
                    jj = floor(posInfo.phi(ii) / 90);
                    kk = find(theta_list_original == posInfo.theta(ii));
                    phi = [linspace(90*jj,90*(jj+1),resNum)'; linspace(90*(jj+1),90*jj,resNum)'];
                    theta = [linspace(theta_list(kk),theta_list(kk),resNum)'; linspace(theta_list(kk+1),theta_list(kk+1),resNum)'];
                    [x,y,z] = sph2cart(phi.*degree, theta.*degree, ones(size(phi)));
                    patchPoints(:,:,ii) = [x,y,z];
            end
        end
    end
    plotPara.x = squeeze(patchPoints(:,1,:));
    plotPara.y = squeeze(patchPoints(:,2,:));
    plotPara.z = squeeze(patchPoints(:,3,:));
end


% create original drp library with inputing angular resolution
function drpLib = createDRPLib(posInfo, ang_res, options)
% this function requires MTEX package for 
    arguments
        posInfo (1,1) struct
        ang_res (1,1) double 
        options.verbose (1,1) logical = 1
        options.crystalSymmetry (1,:) string = 'cubic'
        options.faceting (1,3) double = [1, 0, 0]
        options.fitting_para (1,6) double = [4,1,12,3,3,6]
    end
    cs = crystalSymmetry(options.crystalSymmetry);
    ori = equispacedSO3Grid(cs,'resolution',ang_res);
    nn = length(ori.phi1);
    drpLib.drpList = zeros(nn,length(posInfo.phi));
    drpLib.eulerList = zeros(nn,3);
    for ii = 1:nn
        eulerAngle = [ori(ii).phi1, ori(ii).Phi, ori(ii).phi2]./degree;
        drpLib.drpList(ii,:) = drpSimCone(posInfo, eulerAngle, options.faceting, options.fitting_para);
        drpLib.eulerList(ii,:) = eulerAngle;  % in degree
        if options.verbose
            workbar(ii/nn,sprintf("processing %d / %d DRPs",[ii nn]));
        end
    end
end


% simulate a drp from Euler angles to a vector form
function drpList = drpSimCone(posInfo, eulerAngle, faceting, fitting_para)
    arguments
        posInfo (1,1) struct  % column 1 for phi, column 2 for theta
        eulerAngle (1,3) double = [0,0,0]
        faceting (1,3) double = [1,0,0]
        fitting_para (1,6) double = [1,0.7,25,4,0.8,8]
    end
    
    eu1 = eulerAngle(1);
    eu2 = eulerAngle(2);
    eu3 = eulerAngle(3);
    rot_facet = normr(rotate_facet(eu1,eu2,eu3,faceting));
    i_Main = fitting_para(1);
    i_facet = fitting_para(2);
    sd_Main = fitting_para(3);
    sd_facet = fitting_para(4);

    cauchy = @(p,x) p(1) ./ ((1+((x)./p(2)).^2));
    
    vec_DRP = zeros(3,length(posInfo.phi));
    nn = length(posInfo.phi);
    for idx = 1:nn
        tmp_vec = thph2vec(45+posInfo.theta(idx)/2, posInfo.phi(idx));
        vec_DRP(:,idx) = normr(tmp_vec);  % vec_DRP in size of nx3
    end
    
    drpList = zeros(1,nn);
    % major reflectance peak simulation
    for ii = 1:length(rot_facet)
        ref_a1 = rot_facet(ii,:);
        ref = [0 0 -1] - 2 * dot([0 0 -1], ref_a1) * ref_a1;
        tmp_thph = vec2thph(ref);
        tmp_theta = tmp_thph(1);
        dPh=abs(posInfo.phi - tmp_thph(2));
        dPh=abs(dPh - 360 * (dPh>180)); 
        peakDist=acosd(sind(posInfo.theta)*sind(tmp_theta)+cosd(posInfo.theta)*cosd(tmp_theta).*cosd(dPh));
        drpList = max(drpList, cauchy([i_Main, sd_Main], peakDist));
    end
    
    % great circle band simulation
    for ii = 1:size(rot_facet,1)
        for jj = 1:size(rot_facet,1)
            if ii == jj
                continue
            end
            vec_1 = rot_facet(ii,:);
            vec_2 = rot_facet(jj,:);
            if all(vec_1 - vec_2 < 1e-3) || all(vec_1 + vec_2 < 1e-3)
                continue
            end
            gcnorm = normr(cross(vec_1, vec_2));
            peakDistb=zeros(size(drpList));
            peakDista=zeros(size(drpList));
            bandDist=zeros(size(drpList));
            for mm = 1:nn
                peakDista(mm) = acosd(dot(vec_1, vec_DRP(:,mm)));
                peakDistb(mm) = acosd(dot(vec_2, vec_DRP(:,mm)));
                bandDist(mm) = asind(dot(gcnorm, vec_DRP(:,mm)));
            end
    %         peakDist = bandDist; % + min(peakDista, peakDistb).^2.5;
            drpList = max(drpList, cauchy([i_facet, sd_facet], bandDist));
        end
    end    
end



function indexResult = IndexEngine_sDRM(dataNorm, drpLib, options)
    arguments
        dataNorm
        drpLib
        options.K (1,1) double = 1
    end
    % treatment on input drps: normalization, 
    % drp_in = normalizeVec(drp_in);
    [Idx, D] = knnsearch(drpLib.drpList, dataNorm.drplist, K=options.K);
    indexResult.Euler = drpLib.eulerList(Idx(:,1),:);
    indexResult.Idx = Idx;
    indexResult.distance = D;
    indexResult.eulerMap = reshape(indexResult.Euler,dataNorm.num_x,dataNorm.num_y,3);
    indexResult.idxMap = reshape(indexResult.Idx,dataNorm.num_x,dataNorm.num_y,options.K);
    indexResult.distanceMap = reshape(indexResult.distance(:,1),dataNorm.num_x,dataNorm.num_y);
    fprintf("Orientation indexing finished!\n")
end


% function to show DRP in 3d northern hemisphere
% function plotDRP(drplist, posInfo, options)
%     arguments
%         drplist
%         posInfo
%         options.cMap (1,1) string = "Parula"
%         options.type (1,1) string = "3d"
%         options.caxis (1,2) double = [0,1]
%     end
%     nn = length(posInfo.phi);
%     theta_list_original = sort(unique(posInfo.theta),'ascend');
%     theta_list = zeros(1,length(theta_list_original)+1);
%     for ii = 1:length(theta_list_original)
%         if theta_list_original(ii) <= 70
%             theta_list(ii) = theta_list_original(ii);
%         else
%             theta_list(ii) = (theta_list_original(ii) + theta_list_original(ii-1)) / 2;
%         end
%     end
%     theta_list(end) = 90;
%     hold on
%     for ii = 1:nn
%         % datapoints in the lower range
%         if posInfo.theta(ii) <= 70
%             phi_temp = [posInfo.phi(ii)-5; posInfo.phi(ii)+5; posInfo.phi(ii)+5; posInfo.phi(ii)-5];
%             kk = find(theta_list_original == posInfo.theta(ii));
%             theta_temp = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
%             [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
%             patch(x,y,z,drplist(ii),"EdgeColor","none")
%         else
%             switch find(theta_list_original == posInfo.theta(ii))
%                 case 17
%                     phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(16))),'ascend');
%                     jj = floor(posInfo.phi(ii) / 22.5);
%                     phi_lower = phi_list_temp(phi_list_temp > jj*22.5 & phi_list_temp < (jj+1)*22.5);
%                     phi_temp = [22.5*jj; phi_lower'; 22.5*(jj+1); 22.5*(jj+1); 22.5*jj];
%                     kk = find(theta_list_original == posInfo.theta(ii));
%                     theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
%                     [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
%                     patch(x,y,z,drplist(ii),"EdgeColor","none")
%                 case 18
%                     phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(17))),'ascend');
%                     jj = floor(posInfo.phi(ii) / 30);
%                     phi_lower = phi_list_temp(phi_list_temp > jj*30 & phi_list_temp < (jj+1)*30);
%                     phi_temp = [30*jj; phi_lower'; 30*(jj+1); 30*(jj+1); 30*jj];
%                     kk = find(theta_list_original == posInfo.theta(ii));
%                     theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
%                     [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
%                     patch(x,y,z,drplist(ii),"EdgeColor","none")
%                 case 19
%                     phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(18))),'ascend');
%                     jj = floor(posInfo.phi(ii) / 90);
%                     phi_lower = phi_list_temp(phi_list_temp > jj*90 & phi_list_temp < (jj+1)*90);
%                     phi_temp = [90*jj; phi_lower'; 90*(jj+1); 90*(jj+1); 90*jj];
%                     kk = find(theta_list_original == posInfo.theta(ii));
%                     theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
%                     [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
%                     patch(x,y,z,drplist(ii),"EdgeColor","none")
%             end
%         end
%     end
%     % axis and plot setting
%     axis equal
%     colormap(gca,options.cMap)
%     warning('off')
%     set(gca,"Position",[0 0 1 1],"Visible","off")
%     clim(gca,options.caxis)
%     warning('on')
% end
% 

function plotDRP(drp,plotPara,options)
    arguments
        drp
        plotPara
        options.cMap (1,1) string="jet"
        options.format (1,1) string="3d"
        options.caxis (1,2) double=[0 1]
    end
    x3 = plotPara.x;
    y3 = plotPara.y;
    z3 = plotPara.z;
    if options.format == "3d"
        patch(x3,y3,z3,drp,"EdgeColor","none")
    elseif options.format == "2d"
        x2 = x3./(1+z3);
        y2 = y3./(1+z3);
        patch(x2,y2,drp,"EdgeColor","none")
    end
    axis equal
    colormap(gca,options.cMap)
    warning('off')
    set(gca,"Position",[0 0 1 1],"Visible","off")
    clim(gca,options.caxis)
    warning('on')

end

% apply kernel smoothing on DRM original dataset
function datakernel = kernelSmooth(dataNorm,options)
    arguments
        dataNorm struct
        options.kernelSize (1,1) double = 1
    end
    datakernel = dataNorm;
    % seeking the adjacent neighbors within the preset kernel size
    num_x = dataNorm.num_x;
    num_y = dataNorm.num_y;
    num_datapoints = size(dataNorm.drplist,2);
    dataMap = reshape(dataNorm.drplist, num_x, num_y, num_datapoints);
    switch options.kernelSize
        case 1
            kernel = [0 1 0; 1 1 1; 0 1 0];
        case 2
            kernel = ones(3);
        otherwise
            kernel = 1;
    end
    kernel = kernel / sum(kernel,"all");  % rescale the kernel value
    for ii = 1:num_datapoints
        dataMap(:,:,ii) = conv2(dataMap(:,:,ii), kernel, "same");
    end
    datakernel.drplist = reshape(dataMap,[],num_datapoints);
    fprintf("Kernel processing finished!\n")
end


% select position to show corresponding DRPs
function drp_selected = showSampleDRP(dataNorm, plotPara, options)
    arguments
        dataNorm struct
        plotPara struct
        options.cMap (1,1) string = "jet"
    end
    figure('Name','demo_fig');
    fig_temp = median(dataNorm.drpMap,3);
    imshow(fig_temp,[min(fig_temp,[],"all"),max(fig_temp,[],"all")],'Border','tight');
    colormap("parula")
    [x,y] = ginput;
    % press 'enter' to stop
    nn = length(y);
    y = fix(y);
    x = fix(x);
    close(findobj('type','figure','name','demo_fig'));
    
    figure('Position',[200,200,200*(nn+1),200])
    tiledlayout(1,nn+1,'TileSpacing','tight','Padding','compact')
    nexttile(1)
    imshow(fig_temp,[min(fig_temp,[],"all"),max(fig_temp,[],"all")],'Border','tight')
    hold on
    scatter(x,y,72,'x','k')
    for ii = 1:nn
        text(x(ii)+5,y(ii)+5,int2str(ii),'FontSize',14)
    end
    hold off
    
    drp_selected = zeros(nn,dataNorm.num_pixel);
    
    for ii = 1:nn
        % DRP from measurement 
        nexttile(ii+1)
        x_pos = y(ii);
        y_pos = x(ii);
        drp_selected(ii,:) = squeeze(dataNorm.drpMap(x_pos,y_pos,:));
        plotDRP(drp_selected(ii,:), plotPara,cMap=options.cMap)
    end
end


% function to compare DRP measurement and indexed outcome
function checkDetail = checkIndexResult(dataNorm, indexResult, options)
    arguments
        dataNorm
        indexResult
        options.cMap (1,1) string = "jet"
        options.faceting (1,3) double = [1 0 0]
        options.fitting_para (1,6) double = [1,0.7,25,4,0.8,8]
    end
    indexEulerMap = reshape(indexResult.Euler,dataNorm.num_x,dataNorm.num_y,3);
    figure("Name","Index Result");
    imshow(plot_ipf_map(indexResult.eulerMap),"Border","tight")
    [x,y] = ginput;
    % press 'enter' to stop
    nn = length(y);
    y = fix(y);
    x = fix(x);
    close(findobj("type","figure","name","Index Result"));
    posList = [y, x];

    drp_selected = zeros(nn,dataNorm.num_pixel);
    figure('Position',[100,100,200*(nn),200*2])
    tiledlayout(2,nn,'TileSpacing','tight','Padding','compact')
    for ii = 1:nn
        % DRP from measurement 
        ax = nexttile(ii);
        x_pos = y(ii);
        y_pos = x(ii);
        drp_selected(ii,:) = squeeze(dataNorm.drpMap(x_pos,y_pos,:));
        plotDRP(drp_selected(ii,:), dataNorm.plotPara)
        colormap(ax, options.cMap)
        ax = nexttile(ii+nn);
        plotDRP(drpSimCone(dataNorm.posInfo, squeeze(indexEulerMap(x_pos,y_pos,:)), ...
            options.faceting, options.fitting_para), dataNorm.plotPara);
        colormap(ax, options.cMap)
    end
    
    figure,
    imshow(plot_ipf_map(indexResult.eulerMap),"Border","tight")
    hold on
    scatter(x,y,72,'filled','o','black')
    for ii = 1:nn
        text(x(ii)+5,y(ii)+5,int2str(ii),'FontSize',14,'FontWeight','bold')
    end
    hold off

    checkDetail.posList = posList;
    checkDetail.drpSelected = drp_selected;

end


% plot DRPs based on EBSD ground truth by EBSD Id
function plotDRPfromEBSD(posInfo,ebsd,idList,plotPara,options)
    arguments
        posInfo
        ebsd
        idList (1,:) double
        plotPara (1,1) struct
        options.cMap (1,1) string = "jet"
        options.faceting (1,3) double = [1 0 0]
        options.fitting_para (1,6) double = [1,0.7,25,4,0.8,8]
    end
    nn = length(idList);
    figure('Position',[100,100,200*(nn),200])
    tiledlayout(1,nn,'TileSpacing','tight','Padding','compact')
    for ii = 1:length(idList)
        ax = nexttile(ii);
        rotTemp = ebsd(idList(ii)).rotations;
        eulerTemp = [rotTemp.phi1, rotTemp.Phi, rotTemp.phi2] ./ degree;
        plotDRP(drpSimCone(posInfo, eulerTemp, ...
            options.faceting, options.fitting_para), plotPara, cMap=options.cMap);
    end
end



