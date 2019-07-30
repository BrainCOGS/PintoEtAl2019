function Pinto2019_fig1_wholeTrialOpto(analysisFilePath,summaryFile)

% Pinto2019_fig1_wholeTrialOpto(analysisFilePath,summaryFile)
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

%% Analysis configuration
cfg.perfTh              = 0.6;
cfg.nBins               = 4;
cfg.logRegBins          = linspace(10,200,cfg.nBins+1);
cfg.minNumTrials        = 100;
cfg.minNumTrialsPerf    = 100;
cfg.minNumTrialsLogReg  = 0;
cfg.biasEgCl            = [0 0 0];%analysisParams.areaCl(5,:);
cfg.viewAngEgID         = [25 2];
cfg.biasEgID            = [27 3];
cfg.clustWhat           = {'percCorrect','bias_abs','speed','excessTravel','unusualEvents'};
cfg.clustMaxPC          = 3; % number of PCs to use for clustering
cfg.clustMaxK           = 5; % maximal number of clusters to test
cfg.fixedK              = false; % fix num. clusters?
cfg.distMeasure         = 'eucl'; % 'eucl' or 'correlation'
cfg.clustCl             = [60 179 113; 145 203 55; 10 35 140; 150 110 70; 230 186 0]./255; %; 90 115 210 # 6
cfg.clustOrder          = [3 1 2];
cfg.visGuideExcludeID   = nan;
cfg.perfThTowers        = 0.6;
cfg.perfThCtrk          = 0.8;
cfg.towersCl            = widefieldParams.darkgray; %[0 0 0];
cfg.ctrlCl              = widefieldParams.darkgreen; %[.6 .6 .6];
cfg.memCl               = [149 7 232]./255;
cfg.posBins             = 0:5:300;
cfg.niter               = 1000;
cfg.WFmice              = [40 41 43 45 47 48];

%% get combined log for baseline behavior (opto & wf)
load([analysisFilePath 'owf_concatLog'])

%% get combined log for baseline behavior (opto & wf)
lg.firstTrialofBlock = [true diff(lg.meanPerfBlock)~=0];
lgdirty              = lg;
[lg, info]           = cleanupConcatLog(lg, cfg.minNumTrials);
lg.mouseNames        = mice;
info.mouseNames      = mice;
filters              = [];
files                = {};

%%
lg_visGuide.firstTrialofBlock = [true diff(lg_visGuide.meanPerfBlock)~=0];
lgVisGuidedirty               = lg_visGuide;
lg_visGuide                   = cleanupConcatLog(lg_visGuide, cfg.minNumTrials, 0, 0, 0, 0, .8);
lg_visGuide.mouseNames        = mice;

%%
lg_memGuide.firstTrialofBlock = [true diff(lg_memGuide.meanPerfBlock)~=0];
lgMemGuidedirty               = lg_memGuide;
lg_memGuide                   = cleanupConcatLog(lg_memGuide, cfg.minNumTrials);
lg_memGuide.mouseNames        = mice;

%% visual guide: performance with and without distractors
fieldls = fields(lg_visGuide);
hasDistractors   = (lg_visGuide.trialType == analysisParams.leftCode  & cellfun(@(x)(~isempty(x)),lg_visGuide.cuePos_R)) | ...
                   (lg_visGuide.trialType == analysisParams.rightCode & cellfun(@(x)(~isempty(x)),lg_visGuide.cuePos_L));
for iF = 1:numel(fieldls)
  if strcmpi(fieldls{iF},'mouseNames'); continue; end
  lg_visGuide_noDistractors.(fieldls{iF}) = lg_visGuide.(fieldls{iF})(~hasDistractors & lg_visGuide.currMaze == 12);
  lg_visGuide_distractors.(fieldls{iF})   = lg_visGuide.(fieldls{iF})(hasDistractors  & lg_visGuide.currMaze == 12);
end
lg_visGuide_noDistractors = cleanupConcatLog(lg_visGuide_noDistractors, cfg.minNumTrials, 0, 0, 0, 0, .8);
lg_visGuide_distractors   = cleanupConcatLog(lg_visGuide_distractors, cfg.minNumTrials, 0, 0, 0, 0, .8);

%% cleanup
fieldls = fields(lg_visGuide);
idx     = lg_visGuide.mouseID == cfg.visGuideExcludeID;
for iF = 1:numel(fieldls)
  if strcmpi(fieldls{iF},'mouseNames'); continue; end
  lg_visGuide.(fieldls{iF})(idx) = [];
end

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, info, filters, files);
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

%% view angle (decoding and variance)
fprintf('view angle...\n')

taskComp.towers.viewAngSD     = xMouseViewAngSD(lg,cfg.minNumTrials,cfg.posBins);
taskComp.visGuide.viewAngSD   = xMouseViewAngSD(lg_visGuide,cfg.minNumTrials,cfg.posBins);
taskComp.memGuide.viewAngSD   = xMouseViewAngSD(lg_memGuide,cfg.minNumTrials,cfg.posBins);

%% overall performance and motor indicators for each task (by mouse)
taskComp.towers.percCorrect        = xMousePercCorrect(lg,cfg.minNumTrialsPerf,false);
taskComp.visGuide.percCorrect      = xMousePercCorrect(lg_visGuide,cfg.minNumTrialsPerf,false);
taskComp.visGuide.percCorrect_noDistractors    = xMousePercCorrect(lg_visGuide_noDistractors,cfg.minNumTrialsPerf,false);
taskComp.visGuide.percCorrect_distractors      = xMousePercCorrect(lg_visGuide_distractors,cfg.minNumTrialsPerf,false);
taskComp.memGuide.percCorrect      = xMousePercCorrect(lg_memGuide,cfg.minNumTrialsPerf,false);

taskComp.towers.motorErrs.rate     = xMouseMotorErrors(lgdirty);
taskComp.visGuide.motorErrs.rate   = xMouseMotorErrors(lgVisGuidedirty);
taskComp.memGuide.motorErrs.rate   = xMouseMotorErrors(lgMemGuidedirty);

%% stats
taskComp.towers.motorErrs.mouseAvg   = mean(taskComp.towers.motorErrs.rate);
taskComp.towers.motorErrs.mouseSem   = std(taskComp.towers.motorErrs.rate)./sqrt(numel(mice)-1);
taskComp.visGuide.motorErrs.mouseAvg = mean(taskComp.visGuide.motorErrs.rate);
taskComp.visGuide.motorErrs.mouseSem = std(taskComp.visGuide.motorErrs.rate)./sqrt(numel(mice)-1);
taskComp.memGuide.motorErrs.mouseAvg = mean(taskComp.memGuide.motorErrs.rate);
taskComp.memGuide.motorErrs.mouseSem = std(taskComp.memGuide.motorErrs.rate)./sqrt(numel(mice)-1);

% ANOVA
me   = [taskComp.towers.motorErrs.rate; taskComp.visGuide.motorErrs.rate; taskComp.memGuide.motorErrs.rate];
pc   = [taskComp.towers.percCorrect.mousePC'; taskComp.visGuide.percCorrect.mousePC'; taskComp.memGuide.percCorrect.mousePC'];
va   = [taskComp.towers.viewAngSD.viewAngSD'; taskComp.visGuide.viewAngSD.viewAngSD'; taskComp.memGuide.viewAngSD.viewAngSD'];
idme = [ones(size(taskComp.towers.motorErrs.rate)); 1+ones(size(taskComp.visGuide.motorErrs.rate)); 2+ones(size(taskComp.memGuide.motorErrs.rate))];
idpc = [ones(size(taskComp.towers.percCorrect.mousePC')); 1+ones(size(taskComp.visGuide.percCorrect.mousePC')); 2+ones(size(taskComp.memGuide.percCorrect.mousePC'))];
idva = [ones(size(taskComp.towers.viewAngSD.viewAngSD')); 1+ones(size(taskComp.visGuide.viewAngSD.viewAngSD')); 2+ones(size(taskComp.memGuide.viewAngSD.viewAngSD'))];

[taskComp.stats.p_motorErr,taskComp.stats.anovaTable_motorErr,taskComp.stats.anovaStats_motorErr] = ...
     anovan(me,idme,'varnames','task','display','off');
taskComp.stats.multComp_motorErr = multcompare(taskComp.stats.anovaStats_motorErr,'display','off');

[taskComp.stats.p_percCorrect,taskComp.stats.anovaTable_percCorrect,taskComp.stats.anovaStats_percCorrect] = ...
     anovan(pc,idpc,'varnames','task','display','off');
taskComp.stats.multComp_percCorrect = multcompare(taskComp.stats.anovaStats_percCorrect,'display','off');

[taskComp.stats.p_viewAngSD,taskComp.stats.anovaTable_p_viewAngSD,taskComp.stats.anovaStats_p_viewAngSD] = ...
     anovan(va,idva,'varnames','task','display','off');
taskComp.stats.multComp_p_viewAngSD = multcompare(taskComp.stats.anovaStats_p_viewAngSD,'display','off');

%% Psychometric for each task
fprintf('psychometrics...\n')
taskComp.towers.psych     = psychometricFit(lg.choice,lg.nCues_RminusL,true);
taskComp.visGuide.psych   = psychometricFit(lg_visGuide.choice,lg_visGuide.nCues_RminusL,true);
taskComp.towers.psychPercCorrect     = psychometricFitPercCorrect(lg.choice,lg.trialType,lg.nCues_RminusL,true,0:4:16);
taskComp.visGuide.psychPercCorrect   = psychometricFitPercCorrect(lg_visGuide.choice,lg_visGuide.trialType,lg_visGuide.nCues_RminusL,true,0:4:16);

%% load opto data
thisdir = pwd; 
cd(analysisFilePath)
load fullGrid_lastMaze_whole_wt_perfTh60_nBinLogReg4 lsrPerf
lsrPerf_wt   = lsrPerf;
load fullGrid_visualGuideMaze_whole_vgat_perfTh80_nBinLogReg4 lsrPerf
lsrPerf_t4   = lsrPerf;
load fullGrid_combined_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf 
lsrPerf_comb = lsrPerf;
load fullGrid_memGuide_noTowersMaze_whole_vgat_perfTh60_nBinLogReg4 lsrPerf 
lsrPerf_mem  = lsrPerf;
load fullGrid_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf

%% load dataset diagnosis to save, plus x-expt comparisons
% generated with diagnoseLsrDataset.m
load datasetDiagnosis diagnosis
lsrDiagnosis = diagnosis;

load fullGrid_wholeTrial_xExptComp pvals
xExptStats_fullGrid = pvals;
cd(thisdir)


%% Configure figure and panels

layout       = [ 1  1  2  2  2  ... 
               ; 1  1  3  4  5  ... 
               ; 6  10 11 12 13 ...
               ; 7  14 15 8  9 ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;
stats        = [];

%% Plot rig schematics
iPanel          = pnlRigScheme(fig, iPanel, cfg);

%% Plot maze schematics
iPanel          = pnlMzScheme(fig, iPanel, cfg);

%% Plot Overall Performance 
iPanel         = pnlTaskComp(fig, iPanel, taskComp, cfg, 'percCorrect'); 

%% plot view angle SD
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

pctowers    = taskComp.towers.viewAngSD.viewAngSD;
pcvisguide  = taskComp.visGuide.viewAngSD.viewAngSD;
pcmemguide  = taskComp.memGuide.viewAngSD.viewAngSD;

boxplot([pcvisguide';pcmemguide';pctowers'],[ones(size(pcvisguide))'; 1+ones(size(pcmemguide))';2+ones(size(pctowers))'],...
  'colors',[cfg.ctrlCl; cfg.memCl; cfg.towersCl])

yl = get(axs,'ylim');
text(1.5,yl(2)*.95,sprintf('P = %1.2g',taskComp.stats.p_viewAngSD), ...
     'color','k','horizontalAlignment','center','fontsize',8)
   
set(axs,'xtick',1:3,'xticklabel',{'vis.-guided','mem.-guided','accum.-towers'})
rotateXLabels(axs,45)
ylabel('View angle S.D. (deg)')
box off

%% Plot transparent skull and ROIs w/ grid
iPanel          = iPanel + 1;
iPanel          = pnlSkull(fig, iPanel, cfg);

%% Plot normalized \delta % correct for all mazes
[iPanel,stats]  = pnlEffectSizeComp(fig, iPanel, lsrPerf, lsrPerf_t4, lsrPerf_mem, stats, chance, cfg, taskComp);

%% Plot \delta % correct visually guided maze, with schematics
[iPanel,stats]  = pnlVisGuideMzPerf(fig, iPanel, cfg, lsrPerf_t4, stats);

%% Plot \delta % correct main maze, with schematics
[iPanel,stats]  = pnlMainMzPerf(fig, iPanel, cfg, lsrPerf, stats);

%% Plot \delta % correct mem guided maze, with schematics
[iPanel,stats]  = pnlMemMzPerf(fig, iPanel, cfg, lsrPerf_mem, stats);

%% Save figure to disk
stats.xExpt       = xExptStats_fullGrid;

%% save analysis summary file
taskCompStats = stats;
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'lsrDiagnosis', 'taskComp', 'lsrPerf_mem', 'taskCompStats', ...
                     'lsrPerf','lsrPerf_t4','lsrPerf_wt','xExptStats_fullGrid','-v7.3')
  else
    save(summaryFile,'lsrDiagnosis', 'taskComp', 'lsrPerf_mem', 'taskCompStats', ...
                     'lsrPerf','lsrPerf_t4','lsrPerf_wt','xExptStats_fullGrid','-append')
  end
end

versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% rig schematics
function iPanel = pnlRigScheme(fig,iPanel,cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/rigOpto.png']);
imshow(im); axis off; axis image

end

%% ------------------------------------------------------------------------
%% maze schematics
function iPanel = pnlMzScheme(fig,iPanel,cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/taskSchematics_3tasks.png']);
imshow(im); axis off; axis image

end

%% ------------------------------------------------------------------------
%% transparent skull + grid schematics
function iPanel = pnlSkull(fig,iPanel,cfg)

% skull
iPanel   = iPanel + 1;
axs      = fig.panel(iPanel);
im       = imread([cfg.stockFigPath '/transparentSkull.png']);
[nX,~]   = size(im);
PxlPerMM = 27.5;

hold(axs, 'on')
imshow(im); axis off; axis image
plot([30 30+PxlPerMM],[nX-30 nX-30],'w-','linewidth',3)
text(20+PxlPerMM/2,nX-45,'1 mm','horizontalAlignment','center',...
                     'color','w','fontsize',FormatDefaults.legendFontSize)

% grid
iPanel  = iPanel + 1;
axs     = fig.panel(iPanel);
im      = imread([cfg.stockFigPath '/headplateROIs.png']);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
bregma   = [150 167.5];
data     = load([analysisParams.gridpath 'fullGridBilateral.mat'],'grid'); 
grid     = data.grid;

for iLoc = 1:numel(grid)
  for iSide = 1:2
    x = grid{iLoc}(iSide,1)*PxlPerMM + bregma(2);
    y = grid{iLoc}(iSide,2)*PxlPerMM + bregma(1);
    if iLoc == 5
      plot(x,y,'o','color',analysisParams.lsrShade,'markerfacecolor',analysisParams.lsrCl,'markersize',4)
    else
      plot(x,y,'x','color',analysisParams.lsrShade,'linewidth',.5,'markersize',4)
    end
  end
end
end

%% ------------------------------------------------------------------------
%% delta performance, accumulation maze
function [iPanel,stats] = pnlMainMzPerf(fig, iPanel, cfg, lsrPerf, stats)

% maze schematic
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-30,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)

% brain plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
[pvals,alphac]   = retrievePvals(lsrPerf,'percCorrect','boot');
pvals.effectSize = 100 - cell2mat(arrayfun(@(x)([x; x]),stats.normEffectSize.accumul','UniformOutput',false));
[~,cbar]        = plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,...
                                 -(.01/3).*pvals.effectSize,alphac,1,'es_size',axs,true);
pos             = cbar.Position;
cbar.Units      = get(axs,'units');
cbar.Position   = pos;

title('\Delta Perf. (%) Evidence accumulation','fontsize',12)

stats.gridCoord             = pvals.coord;
stats.mainMaze.perf.alpha   = alphac;
stats.mainMaze.perf.pvals   = pvals.p;%(1:2:end);
stats.mainMaze.perf.effSize = pvals.effectSize;%(1:2:end);

end

%% ------------------------------------------------------------------------
%% delta performance, visually guided maze
function [iPanel, stats] = pnlVisGuideMzPerf(fig, iPanel, cfg, lsrPerf, stats)

% maze schematic
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematicVisualGuide.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)

% brain plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
[pvals,alphac]   = retrievePvals(lsrPerf,'percCorrect','boot');
pvals.effectSize = 100 - cell2mat(arrayfun(@(x)([x; x]),stats.normEffectSize.visGuide','UniformOutput',false));
plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,-(.01/3).*pvals.effectSize,alphac,1,'es_size',axs,true);
title(['\Delta ' sprintf('Perf. (%%) No accumulation\n(Visual guide indicates reward location)')],'fontsize',12)

stats.visGuideMaze.perf.alpha   = alphac;
stats.visGuideMaze.perf.pvals   = pvals.p;
stats.visGuideMaze.perf.effSize = pvals.effectSize;

end

%% ------------------------------------------------------------------------
%% delta performance, accumulation maze
function [iPanel,stats] = pnlMemMzPerf(fig, iPanel, cfg, lsrPerf, stats)

% maze schematic
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic_memNoTowers.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-30,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)

% brain plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
[pvals,alphac]   = retrievePvals(lsrPerf,'percCorrect','boot');
pvals.effectSize = 100 - cell2mat(arrayfun(@(x)([x; x]),stats.normEffectSize.memGuide','UniformOutput',false));
plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,-(.01/3).*pvals.effectSize,alphac,1,'es_size',axs,true);

title('\Delta Perf. (%) Evidence accumulation','fontsize',12)

stats.gridCoord             = pvals.coord;
stats.memMaze.perf.alpha   = alphac;
stats.memMaze.perf.pvals   = pvals.p;
stats.memMaze.perf.effSize = pvals.effectSize;

end

%% ------------------------------------------------------------------------
%% comparison of normalized performance drop across tasks
function [iPanel,stats]  = pnlEffectSizeComp(fig, iPanel, lsrPerf, lsrPerf_t4, lsrPerf_mem, stats, chance, cfg, taskComp)

%%
accumul  = 100 - 100.*cell2mat(cellfun(@(x)(abs(x.lsr.percCorrect-x.ctrl.percCorrect)/(x.ctrl.percCorrect - 100*chance.accumul)),...
                    lsrPerf,'uniformoutput',false));
visGuide = 100 - 100.*cell2mat(cellfun(@(x)(abs(x.lsr.percCorrect-x.ctrl.percCorrect)/(x.ctrl.percCorrect - 100*chance.visGuide)),...
                    lsrPerf_t4,'uniformoutput',false));
memGuide = 100 - 100.*cell2mat(cellfun(@(x)(abs(x.lsr.percCorrect-x.ctrl.percCorrect)/(x.ctrl.percCorrect - 100*chance.memGuide)),...
                    lsrPerf_mem,'uniformoutput',false));  

nROI     = numel(lsrPerf);
datavec  = [visGuide' accumul' memGuide'];
roivec   = repmat((1:nROI)',[1 3]);
taskvec  = [ones(nROI,1) ones(nROI,1)+1 ones(nROI,1)+2];

[stats.normPerf_ANOVA_p,stats.normPerf_ANOVA_table,stats.normPerf_ANOVA_stats] ...
         = anovan(datavec(:),{taskvec(:),roivec(:)},'varnames',{'task','location'},'display','off');

stats.normEffectSize.accumul  = accumul;
stats.normEffectSize.visGuide = visGuide;
stats.normEffectSize.memGuide = memGuide;

% do paired one-tailed stats, only visGuide is not normally distributed
[~,stats.normEffectSize.p_accumulVSmemGuide] = ttest(accumul',memGuide','tail','left'); 
stats.normEffectSize.p_accumulVSmemGuide     = signrank(accumul',visGuide','tail','left'); 
stats.normEffectSize.p_accumulVSmemGuide     = signrank(memGuide',visGuide','tail','left'); 

%
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

% figure; hold on; axs = gca;
plot([.65 3.35],[100 100],'--','linewidth',.25,'color',[.7 .7 .7])
plot([.65 3.35],[0 0],'--','linewidth',.25,'color',[.7 .7 .7])
plot(1:3,[visGuide; memGuide; accumul],'-','linewidth',.5,'color',[.5 .5 .5])
errbar(1,mean(visGuide),std(visGuide)./sqrt(nROI-1),cfg.ctrlCl,1.5);
errbar(2,mean(memGuide),std(memGuide)./sqrt(nROI-1),cfg.memCl,1.5);
errbar(3,mean(accumul),std(accumul)./sqrt(nROI-1),cfg.towersCl,1.5);
plot(1,mean(visGuide),'.','color',cfg.ctrlCl,'markersize',15);
plot(2,mean(memGuide),'.','color',cfg.memCl,'markersize',15);
plot(3,mean(accumul),'.','color',cfg.towersCl,'markersize',15);
xlim([.75 3.25])
ylim([-2 102])
set(axs,'xtick',1:3,'xticklabel',{'vis guided','mem guided','accum. towers'})
rotateXLabels(axs,45)
ylabel(sprintf('Norm. laser on perf.\n(%% of laser off)'))


%
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

accumul_pc  = taskComp.towers.percCorrect.mousePC';  
visGuide_pc = taskComp.visGuide.percCorrect.mousePC';  
memGuide_pc = taskComp.memGuide.percCorrect.mousePC';  
pc          = [nanmean(accumul_pc)*ones(size(accumul)) nanmean(visGuide_pc)*ones(size(visGuide)) nanmean(memGuide_pc)*ones(size(memGuide))]';
lsr         = [accumul visGuide memGuide]';


plot(pc,lsr,'+','color',[.7 .7 .7])
errbar(nanmean(accumul_pc),nanmean(accumul),nanstd(accumul)./sqrt(numel(accumul)-1),cfg.towersCl,1.5,0,'none');
errbar(nanmean(accumul_pc),nanmean(accumul),nanstd(accumul_pc)./sqrt(numel(accumul_pc)-1),cfg.towersCl,1.5,1,'none');
errbar(nanmean(visGuide_pc),nanmean(visGuide),nanstd(visGuide)./sqrt(numel(visGuide)-1),cfg.ctrlCl,1.5,0,'none');
errbar(nanmean(visGuide_pc),nanmean(visGuide),nanstd(visGuide_pc)./sqrt(numel(visGuide_pc)-1),cfg.ctrlCl,1.5,1,'none');
errbar(nanmean(memGuide_pc),nanmean(memGuide),nanstd(memGuide)./sqrt(numel(memGuide)-1),cfg.memCl,1.5,0,'none');
errbar(nanmean(memGuide_pc),nanmean(memGuide),nanstd(memGuide_pc)./sqrt(numel(memGuide_pc)-1),cfg.memCl,1.5,1,'none');

thisfit = fit(pc(:),lsr(:),'poly1');
hold on; plot(60:100,feval(thisfit,60:100),'k-');
plot(60:100,predint(thisfit,60:100),'k--')
[r,p]=corr(pc(:),lsr(:));

text(65,95,sprintf('r = %1.2f\np = %1.2g',r,p))
xlim([60 100]); ylim([-3 103])
xlabel('Avg. overall perf. (% correct)')
ylabel('Norm. laser on perf. (%)')


end

%% -------------------------------------------------------------------
%% task comparison
function iPanel = pnlTaskComp(fig, iPanel, taskComp, cfg, plotWhat)

switch plotWhat
  case 'percCorrect'
    towers   = taskComp.towers.percCorrect.mousePC'; 
    visguide = taskComp.visGuide.percCorrect.mousePC';
    memguide = taskComp.memGuide.percCorrect.mousePC';
    p        = taskComp.stats.p_percCorrect;
    ylbl     = 'Perf. (% correct)';
    
  case 'evidModulation'
    towers   = taskComp.evidModulation.towers.PC_highMinusLowDelta; 
    visguide = taskComp.evidModulation.visGuide.PC_highMinusLowDelta;
    p        = taskComp.evidModulation.stats.towersVSvisGuide_p;
    ylbl     = sprintf('Perf. (%s %% correct)\n(high - low evidence)','\Delta');
    
  case 'speed'
    towers   = taskComp.motor.towers.speed; 
    visguide = taskComp.motor.visGuide.speed;
    memguide = taskComp.motor.memGuide.speed;
    p        = taskComp.motor.stats.p_speed;
    ylbl     = 'Speed (cm/s)';
    
  case 'turnOnset'
    towers   = taskComp.motor.towers.turnOnset; 
    visguide = taskComp.motor.visGuide.turnOnset;
    memguide = taskComp.motor.memGuide.turnOnset;
    p        = taskComp.motor.stats.p_turnOnset;
    ylbl     = 'Turn onset (cm)';

  case 'accelFinish'
    towers   = taskComp.motor.towers.accel_last50; 
    visguide = taskComp.motor.visGuide.accel_last50;
    memguide = taskComp.motor.memGuide.accel_last50;
    p        = taskComp.motor.stats.p_accel_last50;
    ylbl     = 'Acceleration (last 50 cm)(cm/s^2)';
end

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

switch plotWhat
  case 'evidModulation'
    boxplot([visguide; towers],[ones(size(visguide)); 1+ones(size(towers))],...
            'colors',[cfg.ctrlCl; cfg.towersCl])
    yl = get(axs,'ylim');
    text(1.5,yl(2)*.95,sprintf('p = %1.2g',p),'color','k','horizontalAlignment','center','fontsize',9)
    set(axs,'xtick',1:3,'xticklabel',{'vis.-guided','accum.-towers'})
    rotateXLabels(axs,45)
    ylabel(ylbl)
    
  otherwise
    boxplot([visguide; memguide ; towers],[ones(size(visguide)); 1+ones(size(memguide)); 2+ones(size(towers))],...
            'colors',[cfg.ctrlCl; cfg.memCl; cfg.towersCl])
    yl = get(axs,'ylim');
    text(2,yl(2)*.95,sprintf('p = %1.2g',p),'color','k','horizontalAlignment','center','fontsize',9)
    set(axs,'xtick',1:3,'xticklabel',{'vis.-guided','mem.-guided','accum.-towers'})
    rotateXLabels(axs,45)
    ylabel(ylbl)
end
box off
end

