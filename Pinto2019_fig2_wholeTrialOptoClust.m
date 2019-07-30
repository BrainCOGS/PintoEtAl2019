function Pinto2019_fig2_wholeTrialOptoClust(analysisFilePath,summaryFile,loadData)

% Pinto2019_fig2_wholeTrialOptoClust(analysisFilePath,summaryFile)
% plots task schematics, performance sumary, whole trial inactivation
% analysisFilePath is path for data analysis files to be loaded, (whole trial vgat, whole trial ctrl etc)
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate all behavioral logs for control trials:
%    concatLogsBehavSummary_owf.m (repo: behaviorAnalysis) (will save owf_concatLog.mat)
% 2) To generate inactivation experiment logs:
%    analyzeConcatLsrLog_batch(exptType). expType is a string that sets the
%    type of data (e.g. whole trial, subtrial etc). Relevant ones for this
%    figure are 'wholeTrial', 'visGuide', 'WTctrls'
%    (will save the fullGrid_*.mat files)(repo: behaviorAnalysis)
% 3) To statistically compare between inactivation experiments:
%    analyzeConcatLsrLog_xExptStats.m (will save
%    fullGrid_wholeTrial_xExptComp.mat)(repo: behaviorAnalysis)
%% ------------------------------------------------------------------------

if nargin < 3; loadData = false; end

%% Analysis configuration
cfg.perfTh              = 0.6;
cfg.nBins               = 4;
cfg.logRegBins          = linspace(10,200,cfg.nBins+1);
cfg.minNumTrials        = 0;
cfg.minNumTrialsLogReg  = 0;
cfg.biasEgCl            = [0 0 0];%analysisParams.areaCl(5,:);
cfg.viewAngEgID         = [25 2];
cfg.biasEgID            = [27 3];
cfg.viewAngEgCl         = analysisParams.areaCl([3 5],:);
cfg.viewAngEgLbl        = {'RSC','aM2'};
cfg.clustWhat           = {'percCorrect','bias_abs','speed','excessTravel','unusualEvents'};
cfg.clustMaxPC          = 3; % number of PCs to use for clustering
cfg.clustMaxK           = 5; % maximal number of clusters to test
cfg.fixedK              = false; % fix num. clusters?
cfg.distMeasure         = 'eucl'; % 'eucl' or 'correlation'
cfg.clustCl             = [60 179 113; 145 203 55; 10 35 140; 150 110 70; 230 186 0]./255; %; 90 115 210 # 6
cfg.clustCl2            = [147 247 89; 19 106 234; 230 186 0]./255; %; 90 115 210 # 6
cfg.clustOrder          = [3 1 2];
cfg.visGuideExcludeID   = nan;%22; % this mouse was not used in towers
cfg.perfThTowers        = 0.6;
cfg.perfThCtrk          = 0.8;
cfg.minNumTrials        = 0;
cfg.minNumTrialsLogReg  = 0;
cfg.towersCl            = widefieldParams.darkgray; %[0 0 0];
cfg.ctrlCl              = widefieldParams.darkgreen; %[.6 .6 .6];
cfg.memCl               = [149 7 232]./255;
cfg.posBins             = 0:5:300;
cfg.niter               = 1000;
cfg.WFmice              = [40 41 43 45 47 48];


%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], []);
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

if loadData
  load(summaryFile,'optoClustStats', 'optoClustStats_memGuide','lsrPerf','lsrPerf_mem')
else
  %% load opto data
  thisdir = pwd; 
  cd(analysisFilePath)
  load fullGrid_memGuide_noTowersMaze_whole_vgat_perfTh60_nBinLogReg4 lsrPerf 
  lsrPerf_mem  = lsrPerf;
  load fullGrid_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf

  %% cluster laser effects
  optoClustStats          = clusterByOptoEffect(lsrPerf,cfg);
  optoClustStats_memGuide = clusterByOptoEffect(lsrPerf_mem,cfg);

  %% save analysis summary file
  if ~isempty(summaryFile)
    if isempty(dir(summaryFile))
      save(summaryFile,'optoClustStats', 'optoClustStats_memGuide','-v7.3')
    else
      save(summaryFile,'optoClustStats', 'optoClustStats_memGuide','-append')
    end
  end
end

%% Configure figure and panels

layout       = [ 1  2  3  4  5    ... 
               ; 6  7  8  8  8    ...
               ; 9  10 11 11 11   ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;
stats        = [];

%% Plot M2 example and map: absolute bias
[iPanel,stats]  = pnlBias(fig, iPanel, cfg, lsrPerf, stats,true,'towers');

%% Plot speed map
[iPanel,stats]  = pnlSpeed(fig, iPanel, cfg, lsrPerf, stats,'towers');

%% Plot M2 example and map: absolute bias, mem guide
[iPanel,stats]  = pnlBias(fig, iPanel, cfg, lsrPerf_mem, stats,false,'memGuide');

%% Plot speed map, mem guide
[iPanel,stats]  = pnlSpeed(fig, iPanel, cfg, lsrPerf_mem, stats,'memGuide');

%% Plot clustering, towers
[iPanel,stats]  = pnlClustOpto(fig, iPanel, cfg, optoClustStats, false, stats);

%% Plot clustering, mem guide
[iPanel,stats]  = pnlClustOpto(fig, iPanel, cfg, optoClustStats_memGuide, true, stats);

%% fig source data
fig2stats = stats;
save(summaryFile,'fig2stats','-append')

%% Save figure to disk
versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end


%% ------------------------------------------------------------------------
%% abs bias, accumulation maze, eg + map
function [iPanel, stats] = pnlBias(fig, iPanel, cfg, lsrPerf, stats,plotEg,whichData)

if plotEg
% example side bias
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

plot([.5 2.5],[0 0],'--','color',FormatDefaults.lightGray)

mc  = lsrPerf{cfg.biasEgID(1)}.ctrl.bias_mouse;
ml1 = lsrPerf{cfg.biasEgID(1)}.lsr.bias_mouse;
ml2 = lsrPerf{cfg.biasEgID(2)}.lsr.bias_mouse;
plot([ones(size(mc)); 2*ones(size(ml1)); 3*ones(size(ml2))],[mc; ml1; ml2],'k-')

stats.panelA_left.xaxis = [ones(size(mc)); 2*ones(size(ml1)); 3*ones(size(ml2))];
stats.panelA_left.yaxis = [mc; ml1; ml2];

   
xlim([.5 2.5]); ylim([-100 100])
data     = load([analysisParams.gridpath 'fullGridBilateral.mat'],'grid'); 
grid     = data.grid;
loc1 = sprintf('ML %1.1f AP%1.1f (aM2)',grid{cfg.biasEgID(1)}(2,1),grid{cfg.biasEgID(1)}(2,2));
loc2 = sprintf('ML %1.1f AP%1.1f (aM2)',grid{cfg.biasEgID(2)}(2,1),grid{cfg.biasEgID(2)}(2,2));
set(axs,'xtick',1:3,'xticklabel',{'lsr off',loc1,loc2})
t = text(.7,15,'Right','color',FormatDefaults.lightGray,'horizontalAlignment','left','fontsize',8);
set(t,'rotation',90);
t = text(.7,-15,'Left','color',FormatDefaults.lightGray,'horizontalAlignment','right','fontsize',8);
set(t,'rotation',90);
rotateXLabels(axs,30);

stats.panelA_left.xlabels = {'lsr off',loc1,loc2};
stats.panelA_left.ylabel  = 'Right side bias (%)';

end
%%
% brain plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
[pvals,alphac]  = retrievePvals(lsrPerf,'bias_abs','mouse');
plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,false);
title('\Delta |bias| (%)','fontsize',12)

switch whichData
  case 'towers'
    stats.mainMaze.bias_abs.alpha   = alphac;
    stats.mainMaze.bias_abs.pvals   = pvals.p(1:2:end);
    stats.mainMaze.bias_abs.effSize = pvals.effectSize(1:2:end);
    
    stats.panelA_right.gridCoord    = pvals.coord;
    stats.panelA_right.effectSize   = pvals.effectSize;
    stats.panelA_right.pvals        = pvals.p;
    stats.panelA_right.dataLabel    = '\Delta |bias| (%), Accum.-towers task';
  case 'memGuide'
    stats.memMaze.bias_abs.alpha   = alphac;
    stats.memMaze.bias_abs.pvals   = pvals.p(1:2:end);
    stats.memMaze.bias_abs.effSize = pvals.effectSize(1:2:end);
    
    stats.panelC.gridCoord    = pvals.coord;
    stats.panelC.effectSize   = pvals.effectSize;
    stats.panelC.pvals        = pvals.p;
    stats.panelC.dataLabel    = '\Delta |bias| (%), Mem.-guided task';
end

end

%% ------------------------------------------------------------------------
%% speed, accumulation maze, eg + map
function [iPanel, stats] = pnlSpeed(fig, iPanel, cfg, lsrPerf, stats, whichData)

% brain plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
[pvals,alphac]  = retrievePvals(lsrPerf,'speed','boot');
plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,false);
title('Speed (% baseline)','fontsize',12)

switch whichData
  case 'towers'
    stats.mainMaze.speed.alpha   = alphac;
    stats.mainMaze.speed.pvals   = pvals.p(1:2:end);
    stats.mainMaze.speed.effSize = pvals.effectSize(1:2:end);
    
    stats.panelB.gridCoord    = pvals.coord;
    stats.panelB.effectSize   = pvals.effectSize;
    stats.panelB.pvals        = pvals.p;
    stats.panelB.dataLabel    = 'Speed (% baseline), Accum.-towers task';
  case 'memGuide'
    stats.memMaze.speed.alpha   = alphac;
    stats.memMaze.speed.pvals   = pvals.p(1:2:end);
    stats.memMaze.speed.effSize = pvals.effectSize(1:2:end);
    
    stats.panelD.gridCoord    = pvals.coord;
    stats.panelD.effectSize   = pvals.effectSize;
    stats.panelD.pvals        = pvals.p;
    stats.panelD.dataLabel    = 'Speed (% baseline), Mem.-guided task';
end

end

%% ------------------------------------------------------------------------
%% clustering of inactivation effects
function [iPanel, figstats] = pnlClustOpto(fig, iPanel, cfg, stats, memFlag, figstats)

%% dendrogram
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);

leafOrder     = optimalleaforder(stats.lnk,stats.distMat); 
if memFlag
  dh          = dendrogram(stats.lnk,'colorthreshold',0.6,'reorder',leafOrder,'labels',{});
  cfg.clustCl = cfg.clustCl2;
  
  figstats.panelH.linkage   = stats.lnk;
  figstats.panelH.leafOrder = leafOrder;
  figstats.panelH.clusterID = stats.clustID;
  figstats.panelH.ylabel    = 'Eucl. distance';
else
  dh          = dendrogram(stats.lnk,'colorthreshold',0.35,'reorder',leafOrder,'labels',{}); % 'orientation','left',

  % reorder cluster numbering for plotting
  kIDold = stats.clustID;
  if ~isempty(cfg.clustOrder)
    for iK = 1:stats.nK
      stats.clustID(kIDold == cfg.clustOrder(iK)) = iK;
    end
  end
  if stats.nK < 4
    cfg.clustCl = cfg.clustCl(1:2:stats.nK*2,:);
  end
  
  figstats.panelE.linkage   = stats.lnk;
  figstats.panelE.leafOrder = leafOrder;
  figstats.panelE.clusterID = stats.clustID;
  figstats.panelE.ylabel    = 'Eucl. distance';
end

% change colors
hold(axs, 'on')
lineColors = cell2mat(get(dh,'Color'));
colorList  = unique(lineColors, 'rows');

for color = 2:size(colorList,1)
  idx                = ismember(lineColors, colorList(color,:), 'rows');
  lineColors(idx, :) = repmat(cfg.clustCl(color-1,:),sum(idx),1);
end
%// Apply the new colors to the chart's line objects (line by line)
for line = 1:numel(dh)
   dh(line).Color     = lineColors(line,:);
   dh(line).LineWidth = 1;
end
set(axs,'xtick',[],'xcolor','w')
ylabel('Eucl. dist.')


%% clusters on brain
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

xoffset  = 60;
ctBrain  = imread(analysisParams.ctBrainPath);
ctBrain  = repmat(ctBrain(:,xoffset+1:size(ctBrain,2)-xoffset,1),[1 1 3]);
ctBrain  = ctBrain+20;
ctBrain(ctBrain>255) = 255;
pxl(:,1) = (stats.ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;
pxl(:,2) = -(stats.AP.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(1);
pxl(:,3) = (-stats.ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;


imshow(ctBrain)
for iK = 1:stats.nK
  plot(pxl(stats.clustID==iK,1),pxl(stats.clustID==iK,2),...
       '.','markersize',10,'color',cfg.clustCl(iK,:))
  plot(pxl(stats.clustID==iK,3),pxl(stats.clustID==iK,2),...
       '.','markersize',10,'color',cfg.clustCl(iK,:))
end
axis off; axis ij; axis image

if memFlag
  figstats.panelI.gridCoordinates_ML   = stats.ML;
  figstats.panelI.gridCoordinates_AP   = stats.AP;
  figstats.panelI.clusterID            = stats.clustID;
else
  figstats.panelF.gridCoordinates_ML   = stats.ML;
  figstats.panelF.gridCoordinates_AP   = stats.AP;
  figstats.panelF.clusterID            = stats.clustID;
end

%% avg effect by cluster
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')
xt = zeros(1,size(stats.dataMat,2));
for iPt = 1:size(stats.dataMat,2)
  xt(iPt) = (iPt-1)*stats.nK + floor(size(stats.dataMat,2)/2) + (iPt-1);
  for iK = 1:stats.nK
    x = (iPt-1)*stats.nK + iK + (iPt-1);
    if iPt == 1
      lh(iK)     = bar(x,100*mean(stats.dataMat(stats.clustID==iK,iPt)),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:));
      lh_lbl{iK} = ['clust ' num2str(iK)];
    else
      bar(x,100*mean(stats.dataMat(stats.clustID==iK,iPt)),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:))
    end
    errbar(x,100*mean(stats.dataMat(stats.clustID==iK,iPt)), ...
               std(100*stats.dataMat(stats.clustID==iK,iPt))./sqrt(sum(stats.clustID==iK)-1),cfg.clustCl(iK,:),.75)
  end
  yl = get(axs,'ylim');
  text(xt(iPt),yl(2)*.98,sprintf('P = %1.2g',stats.ANOVA(iPt).pval)) % print p val
end

legend(lh,lh_lbl,'location','best'); legend('boxoff')
set(axs, 'xtick', xt, 'xticklabel', stats.labels,'ytick',[-30 0 30])
xlim([0 xt(end)+floor(stats.nK/2)+1])
if memFlag; ylim([-60 60]); else; ylim([-30 30]); end
ylabel(axs, 'Effect size (%)')
rotateXLabels(axs,30)


stats.labels{strcmp(stats.labels,'bias abs')} = 'bias_abs';

if memFlag
  figstats.panelJ.clusterID            = stats.clustID;
  for iData = 1:numel(stats.labels)
    figstats.panelJ.(stats.labels{iData})     = stats.dataMat(:,iData);
  end
  figstats.panelJ.xlabels              = stats.labels;
  figstats.panelJ.ylabel               = 'Effect size (%)';
else
  figstats.panelG.clusterID            = stats.clustID;
  for iData = 1:numel(stats.labels)
    figstats.panelG.(stats.labels{iData})     = stats.dataMat(:,iData);
  end
  figstats.panelG.xlabels              = stats.labels;
  figstats.panelG.ylabel               = 'Effect size (%)';
end

end