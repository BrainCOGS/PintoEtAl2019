function Pinto2019_fig4_WFcorr(analysisFilePath,summaryFile)

% Pinto2019_fig4_WFcorr(analysisFilePath,summaryFile)
% plots correlation analysis results
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To analyze correlations for each recording: dffCorr.m (repo: widefieldImaging)
% 2) To summarize experiments:
%    summarizeROIcorr.m (will save ROIcorrSummary.mat, loaded here)(repo: widefieldImaging)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.corrCMap         = gray;
cfg.corrDiffCMap     = red2blue;
cfg.clustWhat        = 'cc_accumul_wholeTrial_correct'; % cluster on this matrix
cfg.clustMaxPC       = 3; % number of PCs to use for clustering
cfg.clustMaxK        = 4; % maximal number of clusters to test
cfg.clustCl          = [60 179 113; 150 110 70; 230 186 0; 10 35 140; 90 115 210]./255; %  
cfg.xrecCorrBins     = -.5:.1:1;
cfg.perfUseClustWith = {'aM2-L','PPC-L'}; % the cluster that contains this area to avoid numerical kID errors 
cfg.leafOrder        = [2 1 4 3 5 6 15 16 8 7 13 14 10 9 12 11];
cfg.corrType         = 'Spearman';

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/ROIcorrSummary.mat'],'corrSumm');

%% compile stats 
stats        = corrSumm.stats;

%% Configure figure and panels
layout       = [ 1:4; 6 5 7 8; 9:12 ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;


%% correlations: spontaneous, visually guided, towers
fig4stats             = [];
[iPanel,fig4stats]    = pnlAvgCorrs(fig, iPanel, corrSumm, cfg, fig4stats);

%% summary of previous panel
[iPanel,stats,fig4stats]        = pnlAvgROICorrByCondition(fig,iPanel,corrSumm,stats,fig4stats);

%% correlations, clustering
[iPanel,stats,order,fig4stats]  = pnlCorrClust(fig, iPanel, corrSumm, cfg, stats,fig4stats);

%% correlations, towers - visually guided (more decrease between clusters)
[iPanel,stats,fig4stats]  = pnlCorrDiffWithClust(fig, iPanel, corrSumm, cfg, stats, order,fig4stats);

%% correlation between performance and within vs across for each clust
[iPanel,stats,fig4stats]  = pnlCorrWithPerf(fig, iPanel, corrSumm, cfg, stats,fig4stats);

%% avg corr by epoch
[iPanel,stats,fig4stats]  = pnlAvgROICorrByEpoch(fig,iPanel,corrSumm,stats,fig4stats);

%% correlations, towers, hard - easy (more decrease between clusters)
[iPanel,stats,fig4stats]  = pnlCorrDiffByDifficulty(fig, iPanel, corrSumm, cfg, stats, order,fig4stats);

%% (strip down a bit for saving)
corrSumm      = rmfield(corrSumm,'mouse');
WFcorrStats   = stats;

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'corrSumm','WFcorrStats','fig4stats','-v7.3')
  else
    save(summaryFile,'corrSumm','WFcorrStats','fig4stats','-append')
  end
end

%% Save figure to disk
versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% avg ROI correlations for different conditions
function [iPanel,figstats] = pnlAvgCorrs(fig,iPanel,corrSumm,cfg,figstats)

plots   = {'spont_overall','visGuide_overall','accumul_overall'};
plotlbl = {'running in dark','control task','towers task'};
datalbl = {'runningInDark','visGuidedTask','accumTowersTask'};
cc      = [];
for iPlot = 1:numel(plots)
  cc(:,:,iPlot) = nanmean(corrSumm.(['cc_' plots{iPlot}]),3);
  figstats.panelA.(['cc_' datalbl{iPlot}]) = cc(:,:,iPlot);
end
nROI    = size(cc,1);
minc    = 0.35; %min(cc(:))

figstats.panelA.ROIlabels = corrSumm.ROIlbl;

iPanel = iPanel + 1;
for iPlot = 1:numel(plots)
  axs  = fig.panel(iPanel);
  
  imagesc(cc(:,:,iPlot),[minc 1]); colormap(axs,cfg.corrCMap)
  set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
  rotateXLabels(gca,90)
  title(plotlbl{iPlot})
  
  if iPlot == numel(plots)
    cbar              = smallcolorbar(axs);
    pos               = cbar.Position;
    cbar.Units        = get(axs,'units');
    cbar.Position     = pos;
    cbar.Label.String = 'r';
  end
  iPanel = iPanel + 1;
end

iPanel = iPanel - 1;
end

%% ------------------------------------------------------------------------
%% avg ROI correlations condition comparison
function [iPanel,stats,figstats] = pnlAvgROICorrByCondition(fig,iPanel,corrSumm,stats,figstats)

plots   = {'spont_overall','visGuide_wholeTrial_correct','accumul_wholeTrial_correct'};
plotlbl = {'running in dark','control task','towers task'};
datalbl = {'runningInDark','visGuidedTask','accumTowersTask'};
cc      = [];
for iPlot = 1:numel(plots)
  cc(:,end+1) = mean(corrSumm.(['conn_avgROIcc_' plots{iPlot}]))';
  figstats.panelB.(['ccAvg_' datalbl{iPlot}]) = cc(:,end);
end
nROI   = size(cc,1);

figstats.panelB.ROIlabels = corrSumm.ROIlbl;

% stats
% one-way ANOVA w/repeated measures + multiple comparisons
stats.corrByCondition.labels = plotlbl;
vecid = [];
for iPlot = 1:numel(plots)
  vecid = [vecid; ones(nROI,1)*iPlot];
end
[stats.corrByCondition.anovaP,stats.corrByCondition.anovaTable,stats.corrByCondition.anovaStats] ...
                               = anovan(cc(:),{vecid,repmat((1:nROI)',[numel(plots) 1])},'display','off');
stats.corrByCondition.multComp = multcompare(stats.corrByCondition.anovaStats,'display','off');

% plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

[areaCl,ROIlbl] = getDefaultROIcl(corrSumm.ROIlbl);


for iROI = 1:numel(ROIlbl)
  plot(1:3,cc(iROI,:),'-','linewidth',.75,'color',areaCl(iROI,:))
end

set(axs, 'xtick', 1:3, 'xticklabel', plotlbl)
xlim(axs,[.75 3.25]); ylim([.4 1])
rotateXLabels(axs, 60)
ylabel('Avg. corr by ROI (r)')

% print stats
plot([1 3],[1   1],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([1 2],[.9 .9],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([2 3],[.8 .8],'-','color',[.3 .3 .3],'linewidth',1.5)
text(2,1.11,sprintf('P = %1.3g',stats.corrByCondition.multComp(2,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(1.5,.96,sprintf('P = %1.3g',stats.corrByCondition.multComp(1,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(2.5,.86,sprintf('P = %1.3g',stats.corrByCondition.multComp(3,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')

end

%% ------------------------------------------------------------------------
%% hiearchical clustering on correlations
function [iPanel,stats,lO,figstats] = pnlCorrClust(fig,iPanel,corrSumm,cfg,stats,figstats)

% select correlation matrix
cc                           = mean(corrSumm.(cfg.clustWhat),3);
stats.clusterCorr.ROIlbl     = corrSumm.ROIlbl;

% pca
[pcload,~,ev]                = pca(cc);
stats.clusterCorr.fracVarPCs = sum(ev(1:cfg.clustMaxPC))/sum(ev);
[stats.clusterCorr.clustID,distMat,lnk,stats.clusterCorr.CHindex] ...
                             = hierClustering(pcload(:,1:cfg.clustMaxPC),'eucl',cfg.clustMaxK,0);
stats.clusterCorr.nK         = numel(unique(stats.clusterCorr.clustID));

if stats.clusterCorr.nK < 4
  cfg.clustCl = cfg.clustCl(1:2:stats.clusterCorr.nK*2,:);
end

% plot dendrogram
iPanel   = iPanel + 1;
[~,lbls] = getDefaultROIcl(stats.clusterCorr.ROIlbl);

th   = 0.7;
axs  = fig.panel(iPanel);
if isempty(cfg.leafOrder)
  lO   = optimalleaforder(lnk,distMat);
else
  lO   = cfg.leafOrder;
end
dh   = dendrogram(lnk,numel(corrSumm.ROIlbl),'labels',lbls,'ColorThreshold', th, 'reorder', lO, 'Orientation', 'right');

figstats.panelC.linkage   = lnk;
figstats.panelC.clusterID = stats.clusterCorr.clustID;
figstats.panelC.leafOrder = lO;

% change colors
hold(axs, 'on')
lineColors = cell2mat(get(dh,'Color'));
colorList  = unique(lineColors, 'rows');

for color = 2:size(colorList,1)
  idx                = ismember(lineColors, colorList(color,:), 'rows');
  lineColors(idx, :) = repmat(cfg.clustCl(color-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:numel(dh)
   dh(line).Color     = lineColors(line,:);
   dh(line).LineWidth = 1;
end
rotateXLabels(gca,90)
% set(axs,'view',[0 180]);

text(1,1.2,sprintf('num clust = %d',stats.clusterCorr.nK))

end

%% ------------------------------------------------------------------------
%% diff between correct trials, accumul vs vis guide, with clusters
function [iPanel,stats,figstats] = pnlCorrDiffWithClust(fig,iPanel,corrSumm,cfg,stats,order,figstats)

%% select correlation matrix, subtract, order by cluster
cc1 = mean(corrSumm.cc_accumul_wholeTrial_correct_nodistractor,3);
cc2 = mean(corrSumm.cc_visGuide_wholeTrial_correct_nodistractor,3);
ccs = cc1 - cc2;
kID = stats.clusterCorr.clustID;

kID         = kID(order);
ccs         = ccs(order,order);
[cl,lbls]   = getDefaultROIcl(stats.clusterCorr.ROIlbl);
lbls        = lbls(order);
cl          = cl(order,:);
nROI        = numel(lbls);

figstats.panelC.delta_cc  = ccs;
figstats.panelC.ROIlabels = lbls;
figstats.panelC.dataLabel = '\Delta r (accum.towers - vis.guided task)';

% diff corr mat
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

imagesc(ccs,[-max(abs(ccs(:))) max(abs(ccs(:)))]); colormap(axs,cfg.corrDiffCMap)
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',lbls,'yticklabel',lbls)
rotateXLabels(gca,90)
title('accumul towers - control task')
axis image
axis tight

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = 'r';

% draw cluster lines
for iK = 1:numel(unique(kID))
  k1 = find(kID == iK, 1, 'first') - .5;
  k2 = find(kID == iK, 1, 'last') + .5;
  plot([k1 k2 k2 k1 k1],[k1 k1 k2 k2 k1],'-','color',cfg.clustCl(iK,:),'linewidth',1.5)
end
axis ij

%% for each ROI, calculate avg corr drop within and across clusters
stats.clusterCorr.deltaCC_withinClust = zeros(nROI,1);
stats.clusterCorr.deltaCC_acrossClust = zeros(nROI,1);

for iROI = 1:nROI
  idx = setdiff(find(kID == kID(iROI)),iROI);
  stats.clusterCorr.deltaCC_withinClust(iROI) = nanmean(ccs(iROI,idx));
  idx = kID ~= kID(iROI);
  stats.clusterCorr.deltaCC_acrossClust(iROI) = nanmean(ccs(iROI,idx));
end

figstats.panelD.delta_cc_withinCluster  = stats.clusterCorr.deltaCC_withinClust;
figstats.panelD.delta_cc_acrossCluster  = stats.clusterCorr.deltaCC_acrossClust;
figstats.panelD.ROIlabels               = lbls;
figstats.panelD.ylabel                  = '\Delta r (accum.towers - vis.guided task)';

% test
within = stats.clusterCorr.deltaCC_withinClust;
across = stats.clusterCorr.deltaCC_acrossClust;
if lillietest(within) || lillietest(across)
  stats.clusterCorr.deltaCC_withinVSxClust_test  = 'signrank';
  stats.clusterCorr.deltaCC_withinVSxClust_p     = signrank(within,across);
  stats.clusterCorr.deltaCC_withinClust_p        = signrank(within);
  stats.clusterCorr.deltaCC_acrossClust_p        = signrank(across);
else
  stats.clusterCorr.deltaCC_withinVSxClust_test  = 'ttest';
  [~,stats.clusterCorr.deltaCC_withinVSxClust_p] = ttest(within,across);
  [~,stats.clusterCorr.deltaCC_withinClust_p]    = ttest(within);
  [~,stats.clusterCorr.deltaCC_acrossClust_p]    = ttest(across);
end


% plot results
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot([.75 2.25],[0 0],'--','color',[.7 .7 .7])
for iROI = 1:nROI; plot([1 2],[within(iROI) across(iROI)],'-','color',cl(iROI,:)); end

xlim(axs, [.75 2.25])
ylim(axs, [-.12 .05])
set(axs,'ytick',-.1:.05:.05)
yl = get(axs, 'ylim');
text(1.5,yl(2)*.95,sprintf('P = %1.2g',stats.clusterCorr.deltaCC_withinVSxClust_p), ...
     'horizontalAlignment','center','color','k','fontsize',9)
set(axs,'xtick',1:2,'xticklabel',{'within clust','across clust'})
ylabel('\Delta Corr (accumul - control task)')
rotateXLabels(gca,60)

end

%% correlation between performance and correlation
function [iPanel,stats,figstats]  = pnlCorrWithPerf(fig, iPanel, corrSumm, cfg, stats,figstats)

%% find cluster number
[~,lbls]  = getDefaultROIcl(stats.clusterCorr.ROIlbl);
idx1      = strcmpi(lbls,cfg.perfUseClustWith{1});
kID1      = stats.clusterCorr.clustID(idx1);
ROIidx1   = stats.clusterCorr.clustID == kID1;
idx2      = strcmpi(lbls,cfg.perfUseClustWith{2});
kID2      = stats.clusterCorr.clustID(idx2);
ROIidx2   = stats.clusterCorr.clustID == kID2;

% get correlation within ROIs in this cluster
nmice    = numel(corrSumm.mouse);
ccvals   = [];
for iMouse = 1:nmice
  nrec   = size(corrSumm.mouse(iMouse).cc_accumul_overall,3);
  for iRec = 1:nrec
    thiscc = corrSumm.mouse(iMouse).cc_accumul_overall(:,:,iRec);
    vals   = thiscc(ROIidx1,ROIidx2);
    ccvals(end+1,:) = mean(vals(:));
  end
end

figstats.panelE.xaxis  = corrSumm.allrecs_perf;
figstats.panelE.yaxis  = ccvals;
figstats.panelE.xlabel = 'Session performance (accum.-towers task, % correct)';
figstats.panelE.ylabel = 'Avg r, frontal vs. parietal ctx (clust 4 vs. clust 2)';

% correlate with performance
perf = corrSumm.allrecs_perf;
[stats.perfVSwithinClust_corr,stats.perfVSwithinClust_p] = corr(perf,ccvals,'type',cfg.corrType);

% plot
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot(perf*100,ccvals,'o','color',cfg.clustCl(4,:))
xlim([50 90])
ylim([.5 .9])
xlabel('Session perf (% correct)')
ylabel(sprintf('Avg. r between parietal and frontal\ncortex (clusts 2 vs 4)'))
text(55,.6,sprintf('r = %1.2f\nP = %1.2g',stats.perfVSwithinClust_corr,stats.perfVSwithinClust_p), ...
     'color','k','fontsize',8)


end

%% ------------------------------------------------------------------------
%% avg ROI correlations condition comparison
function [iPanel,stats,figstats] = pnlAvgROICorrByEpoch(fig,iPanel,corrSumm,stats,figstats)

plots   = {'accumul_cueHalf1_correct','accumul_cueHalf2_correct','accumul_delay_correct'};
plotlbl = {'cue half 1 (0-1m)','cue half 2 (1-2m)','delay (2-3m)'};
cc      = [];
for iPlot = 1:numel(plots)
  cc(:,end+1) = mean(corrSumm.(['conn_avgROIcc_' plots{iPlot}]))';
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
xlim(axs,[.75 3.25]); ylim([.4 .7])
rotateXLabels(axs, 60)
ylabel('Avg. r by ROI')

% print stats
plot([1 3],[.8   .8],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([1 2],[.75 .75],'-','color',[.3 .3 .3],'linewidth',1.5)
plot([2 3],[.7 .7],'-','color',[.3 .3 .3],'linewidth',1.5)
text(2,.83,sprintf('P = %1.3g',stats.corrByEpoch.multComp(2,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(1.5,.78,sprintf('P = %1.3g',stats.corrByEpoch.multComp(1,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')
text(2.5,.73,sprintf('P = %1.3g',stats.corrByEpoch.multComp(3,end)),'color',[.3 .3 .3],'fontsize',8,'horizontalAlignment','center')

figstats.panelF.cue1stHalf  = cc(:,1);
figstats.panelF.cue2ndHalf  = cc(:,2);
figstats.panelF.delay       = cc(:,3);
figstats.panelF.ROIlabels   = ROIlbl;
figstats.panelF.ylabel      = 'Avg r by ROI';

end

%% ------------------------------------------------------------------------
%% diff between correct trials, accumul vs vis guide, with clusters
function [iPanel,stats,figstats] = pnlCorrDiffByDifficulty(fig,iPanel,corrSumm,cfg,stats,order,figstats)

%% select correlation matrix, subtract, order by cluster
cc1 = mean(corrSumm.cc_accumul_wholeTrial_hard,3);
cc2 = mean(corrSumm.cc_accumul_wholeTrial_easy,3);
ccs = cc1 - cc2;
kID = stats.clusterCorr.clustID;

kID         = kID(order);
ccs         = ccs(order,order);
[cl,lbls]   = getDefaultROIcl(stats.clusterCorr.ROIlbl);
lbls        = lbls(order);
cl          = cl(order,:);
nROI        = numel(lbls);

figstats.panelG.delta_cc  = ccs;
figstats.panelG.ROIlabels = lbls;
figstats.panelG.dataLabel = '\Delta r (hard - easy trials, accum.towers task)';

% diff corr mat
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

imagesc(ccs,[-max(abs(ccs(:))) max(abs(ccs(:)))]); colormap(axs,cfg.corrDiffCMap)
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',lbls,'yticklabel',lbls)
rotateXLabels(gca,90)
title('hard - easy trials')
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

%% for each ROI, calculate avg corr drop within and across clusters
stats.clusterCorr.deltaCC_withinClust_difficulty = zeros(nROI,1);
stats.clusterCorr.deltaCC_acrossClust_difficulty = zeros(nROI,1);

for iROI = 1:nROI
  idx = setdiff(find(kID == kID(iROI)),iROI);
  stats.clusterCorr.deltaCC_withinClust_difficulty(iROI) = nanmean(ccs(iROI,idx));
  idx = kID ~= kID(iROI);
  stats.clusterCorr.deltaCC_acrossClust_difficulty(iROI) = nanmean(ccs(iROI,idx));
end

% test
within = stats.clusterCorr.deltaCC_withinClust_difficulty;
across = stats.clusterCorr.deltaCC_acrossClust_difficulty;
if lillietest(within) || lillietest(across)
  stats.clusterCorr.deltaCC_withinVSxClust_test_difficulty  = 'signrank';
  stats.clusterCorr.deltaCC_withinVSxClust_p_difficulty     = signrank(within,across);
  stats.clusterCorr.deltaCC_withinClust_p_difficulty        = signrank(within);
  stats.clusterCorr.deltaCC_acrossClust_p_difficulty        = signrank(across);
else
  stats.clusterCorr.deltaCC_withinVSxClust_test_difficulty  = 'ttest';
  [~,stats.clusterCorr.deltaCC_withinVSxClust_p_difficulty] = ttest(within,across);
  [~,stats.clusterCorr.deltaCC_withinClust_p_difficulty]    = ttest(within);
  [~,stats.clusterCorr.deltaCC_acrossClust_p_difficulty]    = ttest(across);
end


% plot results
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

for iROI = 1:nROI; plot([1 2],[within(iROI) across(iROI)],'-','color',cl(iROI,:)); end

xlim(axs, [.75 2.25])
yl = get(axs, 'ylim');
text(1.5,yl(2)*.95,sprintf('p = %1.2g',stats.clusterCorr.deltaCC_withinVSxClust_p_difficulty), ...
     'horizontalAlignment','center','color','k','fontsize',9)
set(axs,'xtick',1:2,'xticklabel',{'within clust','across clust'})
ylabel('\Delta Corr (hard - easy trials)')
rotateXLabels(gca,60)

figstats.panelH.delta_cc_withinCluster  = within;
figstats.panelH.delta_cc_acrossCluster  = across;
figstats.panelH.ROIlabels               = lbls;
figstats.panelH.ylabel                  = '\Delta r (hard - easy trials, accum.towers task)';

end
