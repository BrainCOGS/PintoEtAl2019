function Pinto2019_figS7_WFcorr(analysisFilePath,summaryFile)

% Pinto2019_figS7_WFcorr(analysisFilePath,summaryFile)
% plots correlation analysis results vs. inter-ROI distance and
% correlation analysis based on inactivation grid
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To analyze correlations for each recording: dffCorr.m (repo: widefieldImaging)
% 2) To summarize experiments:
%    summarizeROIcorr.m (will save ROIcorrSummary.mat, loaded here)(repo: widefieldImaging)
% 3) ROIs previously generated as described in owf_WFdynamics.m
% 4) To analyze inactivation-grid-based correlations for each recording:
%    dffCorr_grid.m (repo: widefieldImaging)
% 5) To summarize grid correlations: 
%    summarizeROIcorr_grid.m (will save ROIcorrSummary_grid.mat, loaded here)(repo: imaging)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.corrCMap         = gray;
cfg.corrDiffCMap     = red2blue;
cfg.clustWhat        = 'cc_accumul_wholeTrial_correct'; % cluster on this matrix
cfg.clustMaxPC       = 4; % number of PCs to use for clustering
cfg.clustMaxK        = 5; % maximal number of clusters to test
cfg.clustCl          = [60 179 113; 145 203 55; 230 186 0; 150 110 70; 10 35 140; 90 115 210]./255;
cfg.corrBins         = .2:.1:1;
cfg.clustOrder       = [4 5 3 2 1];

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/ROIcorrSummary.mat'],'corrSumm');

%% Configure figure and panels
layout       = [1:3; 4:6; 7:9] ;
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% relationship between distance and correlation
% distance here is defined by the average ROI centroid
% calculate separetly then average over hemispheres
load(summaryFile,'regMaps');
ROIr = regMaps.gcamp_avg.ROI(2:2:end);
ROIl = regMaps.gcamp_avg.ROI(1:2:end);
ROIcentr_r = cellfun(@mean,ROIr,'uniformOutput',false);
ROIcentr_l = cellfun(@mean,ROIl,'uniformOutput',false);
nROI       = numel(ROIcentr_r);
dist_r     = nan(nROI,nROI);
dist_l     = dist_r;

for iR = 1:nROI
  for iC = 1:nROI
    if iR < iC
      dist_r(iR,iC) = sqrt((ROIcentr_r{iR}(1) - ROIcentr_r{iC}(1))^2 + (ROIcentr_r{iR}(2) - ROIcentr_r{iC}(2))^2);
      dist_l(iR,iC) = sqrt((ROIcentr_l{iR}(1) - ROIcentr_l{iC}(1))^2 + (ROIcentr_l{iR}(2) - ROIcentr_l{iC}(2))^2);
    end
  end
end

fact = widefieldParams.pxlPerMM / widefieldParams.dsFactor;
dist = mean(cat(3,dist_r,dist_l),3)./fact;
dist = dist(:);

types = {'spont','visGuide','accumul'};
lbls  = {'Spontaneous','Control','Towers'};
for iType = 1:numel(types)
  thiscc = eval(['corrSumm.cc_' types{iType} '_overall']);
  cc_r   = mean(thiscc(2:2:end,2:2:end,:),3);
  cc_l   = mean(thiscc(1:2:end,1:2:end,:),3);
  cc     = mean(cat(3,cc_r,cc_l),3);
  cc     = cc(:);
  [WFcorrStatsVSdist.ccVSdist_r.(types{iType}),WFcorrStatsVSdist.ccVSdist_p.(types{iType})] ...
         = corr(dist(~isnan(dist)),cc(~isnan(dist)));
       
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  hold(axs,'on')
  plot(dist,cc,'o','color',[.5 .5 .5])
  xlim([0 8]); ylim([.35 1])
  xlabel('Pairwise ROI dist. (mm)')
  ylabel('Pairwise ROI CC (r)')
  set(axs,'xtick',0:2:8,'ytick',.4:.2:1)
  txt = sprintf('r = %1.2f\np = %1.2g',WFcorrStatsVSdist.ccVSdist_r.(types{iType}),WFcorrStatsVSdist.ccVSdist_p.(types{iType}));
  text(.5,.5,txt,'color','k','fontsize',8)
  title(lbls{iType})
end

%% avg corr by "epochs" target task
[iPanel,stats]  = pnlAvgROICorrByEpoch(fig,iPanel,corrSumm);
WFcorrStatsVSdist.visGuideCCbyEpoch = stats;
clear stats corrSumm

%% grid-based correlations - clustering and relationship to inactivation
load([analysisFilePath '/ROIcorrSummary_grid.mat'],'corrSumm');

% select correlation matrix
cc                           = mean(corrSumm.(cfg.clustWhat),3);
stats.clusterCorr.ROIlbl     = corrSumm.ROIlbl;

% average within hemispheres
% cc_all   = cc;
cc_left  = cc(1:2:end,1:2:end);
cc_right = cc(2:2:end,2:2:end);
cc       = mean(cat(3,cc_left,cc_right),3);

% pca
[pcload,~,ev]                = pca(cc);
stats.clusterCorr.fracVarPCs = sum(ev(1:cfg.clustMaxPC))/sum(ev);
[stats.clusterCorr.clustID,~,stats.clusterCorr.lnk,stats.clusterCorr.CHindex] ...
                             = hierClustering(pcload(:,1:cfg.clustMaxPC),'eucl',cfg.clustMaxK,0);
stats.clusterCorr.nK         = numel(unique(stats.clusterCorr.clustID));

% reorder cluster numbering for plotting
kIDold = stats.clusterCorr.clustID;
if ~isempty(cfg.clustOrder)
  for iK = 1:stats.clusterCorr.nK
    stats.clusterCorr.clustID(kIDold == cfg.clustOrder(iK)) = iK;
  end
end
if stats.clusterCorr.nK < 4
  cfg.clustCl = cfg.clustCl(1:2:stats.nK*2,:);
end

% load inactivation grid and compile 
data     = load([analysisParams.gridpath 'fullGridBilateral.mat'],'grid'); 
grid     = data.grid;
oldID = stats.clusterCorr.clustID;
stats.clusterCorr.clustID = [];
stats.clusterCorr.ML = [];
stats.clusterCorr.AP = [];
for iLoc = 1:numel(grid)
  stats.clusterCorr.ML(end+1,:) = grid{iLoc}(1,1);
  stats.clusterCorr.AP(end+1,:) = grid{iLoc}(1,2);
  stats.clusterCorr.ML(end+1,:) = grid{iLoc}(2,1);
  stats.clusterCorr.AP(end+1,:) = grid{iLoc}(2,2);
  stats.clusterCorr.clustID(end+1,:) = oldID(iLoc);
  stats.clusterCorr.clustID(end+1,:) = oldID(iLoc);
end

% plot clustering results on brain 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

xoffset  = 60;
ctBrain  = imread(analysisParams.ctBrainPath);
ctBrain  = repmat(ctBrain(:,xoffset+1:size(ctBrain,2)-xoffset,1),[1 1 3]);
ctBrain  = ctBrain+20;
ctBrain(ctBrain>255) = 255;
pxl      = [];
pxl(:,1) = (stats.clusterCorr.ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;
pxl(:,2) = -(stats.clusterCorr.AP.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(1);

imshow(ctBrain); 
hold(axs,'on')
for iC = 1:stats.clusterCorr.nK
  plot(pxl(stats.clusterCorr.clustID==iC,1),pxl(stats.clusterCorr.clustID==iC,2),...
    '.','markersize',20,'color',cfg.clustCl(iC,:))
  hold(axs,'on')
end
axis off; axis ij; axis image

% correlation between activity and inactivation effects
load(summaryFile,'optoClustStats')
inactCorr = optoClustStats.corrMat;
inactCorr(triu(true(size(inactCorr)),0)) = nan;
inactCorr = inactCorr(:);
inactCorr(isnan(inactCorr)) = [];

activCorr = cc;
activCorr(triu(true(size(activCorr)),0)) = nan;
activCorr = activCorr(:);
activCorr(isnan(activCorr)) = [];

% plot scatter and overlaid binned
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

plot(inactCorr,activCorr,'o','color',[.7 .7 .7],'markersize',4,'linewidth',.25)
xval = toBinCenters(cfg.corrBins);
yval = zeros(size(xval));
sem  = yval;
for iBin = 1:numel(cfg.corrBins)-1
  idx = inactCorr > cfg.corrBins(iBin) & inactCorr <= cfg.corrBins(iBin+1);
  yval(iBin) = mean(activCorr(idx));
  sem(iBin)  = std(activCorr(idx))/sqrt(sum(idx)-1);
end
[r,p]=corr(inactCorr,activCorr);
errorbar(xval,yval,sem,'k-','linewidth',.75)
text(.25,.9,sprintf('r = %1.2f\np = %1.1g',r,p))
xlim([cfg.corrBins(1) cfg.corrBins(end)]) 
ylim([cfg.corrBins(1) cfg.corrBins(end)])
set(axs,'xtick',.2:.4:1,'xtick',.2:.4:1)
xlabel('Inactivation effect similarity (r)')
ylabel('Activity similarity (r)')

stats.inactVSactivCorr_r = r;
stats.inactVSactivCorr_p = p;

%% grid results

% sort by cluster
[~,order] = sort(stats.clusterCorr.clustID);

% plot towers - visual guide for grid
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

cc1 = mean(corrSumm.cc_accumul_wholeTrial_correct,3);
cc2 = mean(corrSumm.cc_visGuide_wholeTrial_correct,3);
ccs = cc1 - cc2;
kID = stats.clusterCorr.clustID;

kID         = kID(order);
ccs         = ccs(order,order);

imagesc(ccs,[-max(abs(ccs(:))) max(abs(ccs(:)))]); colormap(axs,cfg.corrDiffCMap)
xlabel('ROI'); ylabel('ROI')
title('accumul towers - control task')
axis image
axis tight

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\Delta r';

% draw cluster lines
for iK = 1:numel(unique(kID))
  k1 = find(kID == iK, 1, 'first') - .5;
  k2 = find(kID == iK, 1, 'last') + .5;
  plot([k1 k2 k2 k1 k1],[k1 k1 k2 k2 k1],'-','color',cfg.clustCl(iK,:),'linewidth',1.5)
end
axis ij

% plot hard - easy for grid
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

cc1 = mean(corrSumm.cc_accumul_wholeTrial_hard,3);
cc2 = mean(corrSumm.cc_accumul_wholeTrial_easy,3);
ccs = cc1 - cc2;
kID = stats.clusterCorr.clustID;

kID         = kID(order);
ccs         = ccs(order,order);

imagesc(ccs,[-max(abs(ccs(:))) max(abs(ccs(:)))]); colormap(axs,cfg.corrDiffCMap)
xlabel('ROI'); ylabel('ROI')
title('accumul towers - control task')
axis image
axis tight

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\Delta r';

% draw cluster lines
for iK = 1:numel(unique(kID))
  k1 = find(kID == iK, 1, 'first') - .5;
  k2 = find(kID == iK, 1, 'last') + .5;
  plot([k1 k2 k2 k1 k1],[k1 k1 k2 k2 k1],'-','color',cfg.clustCl(iK,:),'linewidth',1.5)
end
axis ij

%% save analysis summary file
corrSummGrid    = rmfield(corrSumm,'mouse');
gridStats       = stats;

if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'WFcorrStatsVSdist','corrSummGrid','gridStats','-v7.3')
  else
    save(summaryFile,'WFcorrStatsVSdist','corrSummGrid','gridStats','-append')
  end
end

%% Save figure to disk
versionInfo.stats = WFcorrStatsVSdist;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% avg ROI correlations condition comparison
function [iPanel,stats] = pnlAvgROICorrByEpoch(fig,iPanel,corrSumm,stats)

plots   = {'cc_visGuide_cueHalf1_correct','cc_visGuide_cueHalf2_correct','cc_visGuide_delay_correct'};
plotlbl = {'0-1 m','1-2 m','2-3 m'};
cc      = [];
for iPlot = 1:numel(plots)
  thiscc = mean(corrSumm.(plots{iPlot}),3);
  thiscc(logical(eye(size(thiscc)))) = nan;
	thiscc = nanmean(thiscc,2);
  cc(:,end+1) = nanmean(thiscc,2);
end
nROI   = size(cc,1);

% stats
% one-way ANOVA w/repeated measures + multiple comparisons
stats.corrByCondition.labels = plotlbl;
vecid = [];
for iPlot = 1:numel(plots)
  vecid = [vecid; ones(nROI,1)*iPlot];
end
[stats.corrByEpoch.anovaP,stats.corrByEpoch.anovaTable,stats.corrByEpoch.anovaStats] ...
                               = anovan(cc(:),{vecid,repmat((1:nROI)',[numel(plots) 1])},'display','off');
stats.corrByEpoch.multComp = multcompare(stats.corrByEpoch.anovaStats,'display','off');

% plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

[areaCl,ROIlbl] = getDefaultROIcl(corrSumm.ROIlbl);


for iROI = 1:numel(ROIlbl)
  plot(1:3,cc(iROI,:),'-','linewidth',.75,'color',areaCl(iROI,:))
end

set(axs, 'xtick', 1:3, 'xticklabel', plotlbl)
xlim(axs,[.75 3.25]); ylim([.4 .8])
rotateXLabels(axs, 60)
ylabel('Avg. r by ROI')

% print stats
plot([1 3],[.9   .9],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([1 2],[.85 .85],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([2 3],[.8 .8],'-','color',[.3 .3 .3],'linewidth',1.5)
text(2,.83,sprintf('P = %1.3g',stats.corrByEpoch.multComp(2,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(1.5,.78,sprintf('P = %1.3g',stats.corrByEpoch.multComp(1,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(2.5,.73,sprintf('P = %1.3g',stats.corrByEpoch.multComp(3,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')

end