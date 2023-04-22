%% scanning DRM data processing script

[dataRaw, dataBG, posInfo] = loadRawData(...
    AlfromUKRaw20230418181452, AlfromUKBG20230418181452, pixel_coord);

dataNorm.drplist = (dataRaw.drplist - dataBG.drp + 40) * dataBG.gain;
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
ang_res = 3;
drpLib = createDRPLib(posInfo, ang_res*degree, faceting=[1,0,0], ...
    fitting_para=[1,0.7,12,4,0.8,8]);
%%
tic
indexResult = IndexEngine_NewDRM(dataNorm.drplist, drpLib.drpList, drpLib.eulerList);
toc
%% plot indexing results
folderpath = "/Users/chenyangzhu/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Research-2023/scanningDRM/dataProcessing";
offset_values = 0:3:69;
errmisOri_stack = zeros(length(offset_values),n1,n2);
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
    errmisOri_stack(ii,:,:) = misOriMap;
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

%%
% eumap = repositionData(indexResult.Euler, dataNorm.x, dataNorm.y);
% figure, imshow(plot_ipf_map(eumap))
figure(Position=[100 100 800 800])
tiledlayout(4,4,"TileSpacing","compact","padding","compact")
rand_idx = randi(size(drpLib.drpList,1),16);
for idx = 1:length(rand_idx)
    nexttile(idx)
    plotDRP(drpLib.drpList(rand_idx(idx),:), posInfo)
end
figure(Position=[100 100 800 800])
tiledlayout(4,4,"TileSpacing","compact","padding","compact")
rand_idx = randi(size(dataNorm.drplist,1),16);
for idx = 1:length(rand_idx)
    nexttile(idx)
    plotDRP(dataNorm.drplist(rand_idx(idx),:), posInfo)
end
%% quick test of the functions
% err_stack = zeros(24,2);
% for ii = 1:24
%     err_stack(ii,1) = mean(errmisOri_stack(ii,:,:),'all',"omitmissing");
%     err_stack(ii,2) = median(errmisOri_stack(ii,:,:),'all','omitmissing');
% end
% figure, plot(offset_values,err_stack(:,1),'LineWidth',2)
% hold on
% plot(offset_values,err_stack(:,2),'LineWidth',2)
% set(gca,'LineWidth',2,'FontSize',14)
% legend("average error","median error")

% Define file names and duration of each frame
gifImgFolder = "C:\Users\86198\OneDrive - Nanyang Technological University\Research-2023\scanningDRM\dataProcessing";
image_prefix = "err_histo_offset_";
fileNames = strcat(image_prefix,sprintf("%02d.tif",offset_values(1)));
durations = 0.5;

% Initialize GIF file
gifFileName = 'err_histo_offset.gif';
for ii = 1:length(offset_values)
    % Read image
    img = imread(fullfile(gifImgFolder,strcat(image_prefix,sprintf("%02d.tif",offset_values(ii)))));
    [ind, map] = rgb2ind(img, 256);
    % Write image to GIF file
    if ii == 1
        % For first image, create new file with overwrite option
        imwrite(ind, map, fullfile(gifImgFolder,gifFileName), 'gif', 'Loopcount', inf, 'DelayTime', durations, 'WriteMode', 'overwrite');
    else
        % For subsequent images, append to existing file
        imwrite(ind, map, fullfile(gifImgFolder,gifFileName), 'gif', 'WriteMode', 'append', 'DelayTime', durations);
    end
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


% function to show DRP in 3d northern hemisphere
function plotDRP(drplist, posInfo, options)
    arguments
        drplist
        posInfo
        options.cMap (1,1) string = "jet"
        options.type (1,1) string = "3d"
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
    hold on
    for ii = 1:nn
        % datapoints in the lower range
        if posInfo.theta(ii) <= 70
            phi_temp = [posInfo.phi(ii)-5; posInfo.phi(ii)+5; posInfo.phi(ii)+5; posInfo.phi(ii)-5];
            kk = find(theta_list_original == posInfo.theta(ii));
            theta_temp = [theta_list(kk); theta_list(kk); theta_list(kk+1); theta_list(kk+1)];
            [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
            patch(x,y,z,drplist(ii),"EdgeColor","none")
        else
            switch find(theta_list_original == posInfo.theta(ii))
                case 17
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(16))),'ascend');
                    jj = floor(posInfo.phi(ii) / 22.5);
                    phi_lower = phi_list_temp(phi_list_temp > jj*22.5 & phi_list_temp < (jj+1)*22.5);
                    phi_temp = [22.5*jj; phi_lower'; 22.5*(jj+1); 22.5*(jj+1); 22.5*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
                case 18
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(17))),'ascend');
                    jj = floor(posInfo.phi(ii) / 30);
                    phi_lower = phi_list_temp(phi_list_temp > jj*30 & phi_list_temp < (jj+1)*30);
                    phi_temp = [30*jj; phi_lower'; 30*(jj+1); 30*(jj+1); 30*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
                case 19
                    phi_list_temp = sort(unique(posInfo.phi(posInfo.theta == theta_list_original(18))),'ascend');
                    jj = floor(posInfo.phi(ii) / 90);
                    phi_lower = phi_list_temp(phi_list_temp > jj*90 & phi_list_temp < (jj+1)*90);
                    phi_temp = [90*jj; phi_lower'; 90*(jj+1); 90*(jj+1); 90*jj];
                    kk = find(theta_list_original == posInfo.theta(ii));
                    theta_temp = [theta_list(kk)*ones(length(phi_lower)+2,1); theta_list(kk+1); theta_list(kk+1)];
                    [x,y,z] = sph2cart(phi_temp.*degree, theta_temp.*degree, ones(size(phi_temp)));
                    patch(x,y,z,drplist(ii),"EdgeColor","none")
            end
        end
    end
    % axis and plot setting
    axis equal
    colormap(options.cMap)
    set(gca,"Visible","off")
end