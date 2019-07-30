function Pinto2019_figS3_wholeTrialOpto_extraInfo(analysisFilePath,summaryFile)

% Pinto2019_figS3_wholeTrialOpto_extraInfo(analysisFilePath)
% plots performance details for whole-trial inactivation
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate inactivation experiment logs:
%    analyzeConcatLsrLog_batch(exptType). expType is a string that sets the
%    type of data (e.g. whole trial, subtrial etc). Relevant ones for this
%    figure are 'wholeTrial', 'visGuide', 'wholeTrialClust'
%    (will save the fullGrid_*.mat files)(repo: behaviorAnalysis)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.perfTh              = 0.6;
cfg.nBins               = 4;
cfg.logRegBins          = linspace(10,200,cfg.nBins+1);
cfg.minNumTrials        = 0;
cfg.minNumTrialsLogReg  = 0;
cfg.biasEgID            = 1;
cfg.biasEgCl            = analysisParams.areaCl(5,:);
cfg.viewAngEgID         = [25 2];
cfg.viewAngEgCl         = analysisParams.areaCl([3 5],:);
cfg.viewAngEgLbl        = {'RSC','aM2'};
cfg.SPplotWhat_vg       = {'bias_abs','speed','excessTravel','unusualEvents'};
cfg.SPplotWhatLbl_vg    = {'|bias| (%)','speed (%)','excess travel (%)','motor errors (%)'};
cfg.SPplotWhat_towers   = {'excessTravel','unusualEvents'};
cfg.SPplotWhatLbl_towers= {'excess travel (%)','motor errors (%)'};
cfg.dataLbl             = {'Clust 1','Clust 2','Clust 3'};
cfg.clustWhat           = {'percCorrect','bias_abs','speed','excessTravel','unusualEvents'};
cfg.clustMaxPC          = 3; % number of PCs to use for clustering
cfg.clustMaxK           = 5; % maximal number of clusters to test
cfg.fixedK              = false; % fix num. clusters?
cfg.distMeasure         = 'eucl'; % 'eucl' or 'correlation'
cfg.clustCl             = [60 179 113; 145 203 55; 10 35 140; 150 110 70; 230 186 0]./255; 
cfg.clustOrder          = [3 1 2];
cfg.clustTestMaxK       = [5 10 5 5];
cfg.clustTestPC         = [2 3 2 3];
cfg.clustTestDist       = {'eucl';'eucl';'correlation';'correlation'};

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

%% load opto data
thisdir = pwd; 
cd(analysisFilePath)

load fullGrid_visualGuideMaze_whole_vgat_perfTh80_nBinLogReg4 lsrPerf
lsrPerf_t4   = lsrPerf;
load fullGrid_memGuide_noTowersMaze_whole_vgat_perfTh60_nBinLogReg4 lsrPerf
lsrPerf_mem  = lsrPerf;
load fullGrid_clusters_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf 
lsrPerf_clust = lsrPerf([1 3 2]);
load fullGrid_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf

cd(thisdir)

if isempty(summaryFile)
  load(summaryFile,'optoClustStats')
else
  optoClustStats = clusterByOptoEffect(lsrPerf,cfg);
end

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'lsrPerf_clust','-v7.3')
  else
    save(summaryFile,'lsrPerf_clust','-append')
  end
end
%% SUPPL FIGURE: combined psychometrics + log Reg, speed, excess travel, motor events
stats             = [];
layout            = [1  3  4  7     ...
                   ; 2  5  6  8     ...
                   ; 9  10 11 21    ...
                   ; 12 13 14 22    ...
                   ; 15 16 17 23    ...
                   ; 18 19 20 24    ...
                    ];
fig               = PaneledFigure(layout, 'smaller');
stats             = figMaps(fig, lsrPerf, lsrPerf_t4, cfg, stats, lsrPerf_clust, lsrPerf_mem);
versionInfo.stats = stats;

fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% extra maps, psychometrics etc
function stats = figMaps(fig, lsrPerf, lsrPerf_t4, cfg, stats, lsrPerf_clust, lsrPerf_mem)

iPanel = 0;

% maze schematic t4
iPanel = iPanel+1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematicVisualGuide.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)

% maze schematic 
iPanel = iPanel+1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-30,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)

% maps
for iP = 1:numel(cfg.SPplotWhat_vg)
  iPanel = iPanel+1;
  axs    = fig.panel(iPanel);
  [pvals,alphac]  = retrievePvals(lsrPerf_t4,cfg.SPplotWhat_vg{iP},'boot');
  [~, cbar]       = plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs);
  title(['\Delta ' cfg.SPplotWhatLbl_vg{iP}],'fontsize',12)
  
  stats.visGuideMaze.(cfg.SPplotWhat_vg{iP}).alpha   = alphac;
  stats.visGuideMaze.(cfg.SPplotWhat_vg{iP}).pvals   = pvals.p(1:2:end);
  stats.visGuideMaze.(cfg.SPplotWhat_vg{iP}).effSize = pvals.effectSize(1:2:end);
  if iPanel == 1
    pos             = cbar.Position;
    cbar.Units      = get(axs,'units');
    cbar.Position   = pos;
  end

end

% maps 
for iP = 1:numel(cfg.SPplotWhat_towers)
  iPanel = iPanel+1;
  axs    = fig.panel(iPanel);
  [pvals,alphac]  = retrievePvals(lsrPerf,cfg.SPplotWhat_towers{iP},'boot');
  plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,iP==1);
  title(['\Delta ' cfg.SPplotWhatLbl_towers{iP}],'fontsize',12)
  
  stats.mainMaze.(cfg.SPplotWhat_towers{iP}).alpha   = alphac;
  stats.mainMaze.(cfg.SPplotWhat_towers{iP}).pvals   = pvals.p(1:2:end);
  stats.mainMaze.(cfg.SPplotWhat_towers{iP}).effSize = pvals.effectSize(1:2:end);
end

% maps 
for iP = 1:numel(cfg.SPplotWhat_towers)
  iPanel = iPanel+1;
  axs    = fig.panel(iPanel);
  [pvals,alphac]  = retrievePvals(lsrPerf_mem,cfg.SPplotWhat_towers{iP},'boot');
  plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,iP==1);
  title(['\Delta ' cfg.SPplotWhatLbl_towers{iP}],'fontsize',12)
  
  stats.memGuideMaze.(cfg.SPplotWhat_towers{iP}).alpha   = alphac;
  stats.memGuideMaze.(cfg.SPplotWhat_towers{iP}).pvals   = pvals.p(1:2:end);
  stats.memGuideMaze.(cfg.SPplotWhat_towers{iP}).effSize = pvals.effectSize(1:2:end);
end

%% combined psychometrics
% iPanel = iPanel + 1;
clustCl = cfg.clustCl([1 3 5],:);
for iArea = 1:numel(lsrPerf_clust)
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  plotPsychometricCurve(lsrPerf_clust{iArea},axs,clustCl(iArea,:),clustCl(iArea,:),false);
  title(cfg.dataLbl{iArea},'fontsize',13)
end

%% combined logistic
for iArea = 1:numel(lsrPerf_clust)
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  plotRevCorr(lsrPerf_clust{iArea},axs,clustCl(iArea,:),true);
  ylim([-.2 .3])
end

%% combined view angle
for iArea = 1:numel(lsrPerf_clust)
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  plotTrajectory(lsrPerf_clust{iArea},'decode',axs,clustCl(iArea,:))
end

%% different clustering params
clustCfg = cfg;
for iK = 1:numel(cfg.clustTestMaxK)
  clustCfg.clustMaxK   = cfg.clustTestMaxK(iK);
  clustCfg.clustMaxPC  = cfg.clustTestPC(iK);
  clustCfg.distMeasure = cfg.clustTestDist{iK};
  
  clustStats = clusterByOptoEffect(lsrPerf,clustCfg);
  
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  hold(axs, 'on')

  xoffset  = 60;
  ctBrain  = imread(analysisParams.ctBrainPath);
  ctBrain  = repmat(ctBrain(:,xoffset+1:size(ctBrain,2)-xoffset,1),[1 1 3]);
  ctBrain  = ctBrain+20;
  ctBrain(ctBrain>255) = 255;
  pxl(:,1) = (clustStats.ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;
  pxl(:,2) = -(clustStats.AP.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(1);
  pxl(:,3) = (-clustStats.ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;


  imshow(ctBrain)
  for iC = 1:clustStats.nK
    plot(pxl(clustStats.clustID==iC,1),pxl(clustStats.clustID==iC,2),...
         '.','markersize',10,'color',cfg.clustCl(iC,:))
    plot(pxl(clustStats.clustID==iC,3),pxl(clustStats.clustID==iC,2),...
         '.','markersize',10,'color',cfg.clustCl(iC,:))
  end
  axis off; axis ij; axis image
  title(sprintf('num. PC = %d, max K = %d\ndist. measure = %s',clustCfg.clustMaxPC,clustCfg.clustMaxK,clustCfg.distMeasure))
end

end

%% ------------------------------------------------------------------------
%% clustering of inactivation effects
function iPanel = pnlClustOpto(fig, iPanel, cfg, stats)

%% dendrogram
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);

leafOrder  = optimalleaforder(stats.lnk,stats.distMat); 
dh         = dendrogram(stats.lnk,'colorthreshold',0.35,'reorder',leafOrder,'labels',{});

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
    errorbar(x,100*mean(stats.dataMat(stats.clustID==iK,iPt)), ...
               std(100*stats.dataMat(stats.clustID==iK,iPt))./sqrt(sum(stats.clustID==iK)-1),...
               '-','color',cfg.clustCl(iK,:))
  end
  yl = get(axs,'ylim');
  text(xt(iPt),yl(2)*.98,sprintf('P = %1.2g',stats.ANOVA(iPt).pval)) % print p val
end

legend(lh,lh_lbl,'location','best'); legend('boxoff')
set(axs, 'xtick', xt, 'xticklabel', stats.labels,'ytick',[-30 0 30])
xlim([0 xt(end)+floor(stats.nK/2)+1])
ylim([-30 30])
ylabel(axs, 'Effect size (%)')
rotateXLabels(axs,30)

end