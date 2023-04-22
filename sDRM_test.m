
%% Seek data distribution of background subtracted data
prcnum = 70;
drpThres = prctile(normValue,prcnum,2);
valueThres = prctile(normValue,prcnum,"all");
figure, histogram(drpThres,BinLimits=[-40 0],BinMethod='integers', ...
    EdgeColor='k',EdgeAlpha=0.5,FaceColor='#0072BD',FaceAlpha=1)
hold on
line([valueThres, valueThres],[0 2e4],"color",'red',"linewidth",4)
set(gca,"LineWidth",2,"FontSize",14)


%% use off-set to percentage number of 50
offset = prctile(normValue,50,'all');
gainVal = 20;

drpVal = normValue - offset;
% remove value less than 0
drpVal(drpVal < 0) = 0;
% apply gain
drpVal = drpVal * gainVal;
drpVal_m = drpVal;
for ii = 1:size(drpVal_m,1)
    drp_temp = drpVal_m(ii,:);
    drp_temp = drp_temp / prctile(drp_temp,95);
    drp_temp(drp_temp > 1) = 1;
    drpVal_m(ii,:) = drp_temp;
    workbar(ii/size(drpVal_m,1))
end

%%
tic
indexResult = IndexEngine_NewDRM(drpVal_m, drpLib.drpList, drpLib.eulerList);
toc

%% functions used in this script
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


