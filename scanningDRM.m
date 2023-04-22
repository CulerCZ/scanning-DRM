%% scanning DRM data processing script

[dataRaw, dataBG, posInfo] = loadRawData(...
    AlfromUKRaw20230418181452, AlfromUKBG20230418181452, pixel_coord);

dataNorm.drplist = (dataRaw.drplist - dataBG.drp + dataBG.offset) * dataBG.gain;
dataNorm.drplist(dataNorm.drplist < 0) = 0;
dataNorm.x = dataRaw.x;
dataNorm.y = dataRaw.y;
dataNorm.num_x = numel(unique(dataNorm.x));
dataNorm.num_y = numel(unique(dataNorm.y));
dataNorm.drplist = dataNorm.drplist / prctile(dataNorm.drplist,95,"all");
dataNorm.drplist(dataNorm.drplist > 1) = 1;

figure, imagesc(unique(dataNorm.x), unique(dataNorm.y),...
    reshape(median(dataNorm.drplist,2),dataNorm.num_x,dataNorm.num_y))

%% create DRP Library and run dictionary indexing
ang_res = 4;
% drpLib = createDRPLib(posInfo, ang_res*degree, faceting=[1,0,0], ...
%     fitting_para=[1,0.7,12,4,0.8,8]);

indexResult = IndexEngine_NewDRM(dataNorm.drplist, drpLib.drpList, drpLib.eulerList);

%% plot indexing results
folderpath = "/Users/chenyangzhu/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Research-2023/scanningDRM/dataProcessing";
offset_values = 0:3:69;
% errmisOri_stack = zeros(length(offset_values),n1,n2);
for ii = 1:length(offset_values)
    dataNorm.drplist = (dataRaw.drplist - dataBG.drp + offset_values(ii)) * 20;
    dataNorm.drplist(dataNorm.drplist < 0) = 0;
    dataNorm.drplist = dataNorm.drplist / prctile(dataNorm.drplist,95,"all");
    dataNorm.drplist(dataNorm.drplist > 1) = 1;
    indexResult = IndexEngine_NewDRM(dataNorm.drplist, drpLib.drpList, drpLib.eulerList);
    eumap = repositionData(indexResult.Euler, dataNorm.x, dataNorm.y);
    
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
    errmisOri_stack(ii+19,:,:) = misOriMap;
    fprintf("Save error map %02d / %02d offset values.\n",[ii length(offset_values)]);
    % plot pixel-wise indexing error histogram
    f2=figure("name",sprintf("offset %d error histogram",offset_values(ii)));
    histogram(misOriMap,61,'BinLimits',[1 62],...
        'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
    set(gca,'LineWidth',1.5,'FontSize',14)
    xlabel('prediction error (deg)')
    xlim([1 62])
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

%% quick test of the functions
err_stack = zeros(24,2);
for ii = 1:24
    err_stack(ii,1) = mean(errmisOri_stack(ii,:,:),'all',"omitmissing");
    err_stack(ii,2) = median(errmisOri_stack(ii,:,:),'all','omitmissing');
end
figure, plot(offset_values,err_stack(:,1),'LineWidth',2)
hold on
plot(offset_values,err_stack(:,2),'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',14)
legend("average error","median error")

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
    posInfo.phi = pixel_coord(1,:);
    posInfo.phi = rem(270-posInfo.phi+360, 360);  % under current settings
    posInfo.theta = pixel_coord(2,:);
end


% quick plot DRP in polar coordinates
function polarPlotDRP(posInfo, drp, options)
    arguments
        posInfo (1,1) struct
        drp double
        options.scatterSize (1,1) double = 50
    end
    x = cosd(posInfo.theta).*cosd(posInfo.phi)./(1+sind(posInfo.theta));
    y = cosd(posInfo.theta).*sind(posInfo.phi)./(1+sind(posInfo.theta));
    scatter(x,y,options.scatterSize,drp,'filled')
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    set(gca,'visible','off')
    
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
        fitting_para (1,6) double = [1,0.7,12,4,0.8,8]
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



function indexResult = IndexEngine_NewDRM(drp_in, drpList, eulerList, options)
    arguments
        drp_in
        drpList
        eulerList
        options.verbose
    end
    % treatment on input drps: normalization, 
    % drp_in = normalizeVec(drp_in);
    [Idx, D] = knnsearch(drpList, drp_in);
    indexResult.Euler = eulerList(Idx,:);
    indexResult.Idx = Idx;
    indexResult.distance = D;
end


% function to re-position the datapoints as pos_x and pos_y
% BUT, just make it use for now...
function dataRepos = repositionData(data, posx, posy)
    numlayer = size(data,2);
    numx = numel(unique(posx));
    numy = numel(unique(posy));
    % idx_x = (posx-min(posx))/((max(posx)-min(posx))/(numx-1))+1;
    % idx_y = (posy-min(posy))/((max(posy)-min(posy))/(numy-1))+1;
    dataRepos = zeros(numx,numy,numlayer);
    for ii = 1:numlayer
        datalayer_temp = reshape(data(:,ii),numx,numy);
        datalayer_temp(:,1:2:end) = flipud(datalayer_temp(:,1:2:end));
        dataRepos(:,:,ii) = datalayer_temp;
    end
end