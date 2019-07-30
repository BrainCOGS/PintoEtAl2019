function Pinto2019_figS3_subTrialOpto(analysisFilePath,summaryFile)

% Pinto2019_figS3_subTrialOpto(analysisFilePath,summaryFile)
% plots task schematics, performance summary, for sub-trial inactivation
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats are saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate inactivation experiment logs:
%    analyzeConcatLsrLog_batch(exptType). expType is a string that sets the
%    type of data (e.g. whole trial, subtrial etc). Relevant one for this
%    figure is 'epochs' (will save the reducedGrid_*.mat files)(repo: behav)
% 2) To statistically compare between inactivation experiments:
%    analyzeConcatLsrLog_xExptStats.m (will save
%    reducedGrid_subTrial_xExptComp.mat)(repo: behaviorAnalysis)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.dataLbl         = {'V1','PPC','RSC','mM2','aM2'};
cfg.egRegionID      = [1 3 4];
cfg.egRegionCl      = analysisParams.areaCl(cfg.egRegionID,:);
cfg.viewAngEgID     = [1 3 4];
cfg.viewAngEgCl     = analysisParams.areaCl(cfg.viewAngEgID,:);
cfg.logRegPlotCtrl  = false;
cfg.viewAngPlotCtrl = false;

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

%% load opto data, compile summary stats
thisdir = pwd; 
cd(analysisFilePath)
load reducedGrid_lastMaze_cueHalf1_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf
lsrPerf_ch1 = lsrPerf;
load reducedGrid_lastMaze_cueHalf2_vgat_perfTh60_6mW_nBinLogReg4  lsrPerf
lsrPerf_ch2 = lsrPerf;
load reducedGrid_lastMaze_mem_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf
lsrPerf_del = lsrPerf;

load reducedGrid_subTrial_xExptComp pvals
xExptStats_reducedGrid = pvals;
cd(thisdir)

subTrialOptoStats       = summaryStats(lsrPerf_ch1, lsrPerf_ch2, lsrPerf_del, cfg);
subTrialOptoStats.xExpt = xExptStats_reducedGrid;

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'lsrPerf_ch1','lsrPerf_ch2','lsrPerf_del','subTrialOptoStats','-v7.3')
  else
    save(summaryFile,'lsrPerf_ch1','lsrPerf_ch2','lsrPerf_del','subTrialOptoStats','-append')
  end
end

%%
lsrPerf_ch1{1}.info.grid{2}(:,2) = [-1.75;-1.75];
lsrPerf_ch2{1}.info.grid{2}(:,2) = [-1.75;-1.75];

layout       = [ 1  2 ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 1;

axs              = fig.panel(iPanel);
[pvals,alphac]   = retrievePvals(lsrPerf_ch1,'percCorrect','boot',.05);
[~,cbar]        = plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,true);
pos             = cbar.Position;
cbar.Units      = get(axs,'units');
cbar.Position   = pos;

iPanel       = 2;

axs              = fig.panel(iPanel);
[pvals,alphac]   = retrievePvals(lsrPerf_ch2,'percCorrect','boot',.05);
plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,true);


%% Save figure to disk
versionInfo.stats = subTrialOptoStats;
fig.export(codeFile, versionInfo, true, true);


end

%% ------------------------------------------------------------------------
%% summary stats
function stats = summaryStats(lsrPerf_ch1, lsrPerf_ch2, lsrPerf_del, cfg)

% compile stats
stats.dataLbl                                    = cfg.dataLbl;

% overall performance
stats.cueHalf1.percCorrect.diffMean              = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch1);
stats.cueHalf1.percCorrect.diffStd               = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_std),lsrPerf_ch1);
stats.cueHalf1.percCorrect.pval                  = cellfun(@(x)(x.stats.bootPerf.p_percCorrect),lsrPerf_ch1);
[~,stats.cueHalf1.percCorrect.alpha]             = FDR(stats.cueHalf1.percCorrect.pval');
stats.cueHalf1.percCorrect.diffMean_vsDelTow     = cellfun(@(x)(x.cueHalf.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch1);
stats.cueHalf1.percCorrect.diffStd_vsDelTow      = cellfun(@(x)(x.cueHalf.stats.bootPerf.diff_percCorrect_std),lsrPerf_ch1);
stats.cueHalf1.percCorrect.pval_vsDelTow         = cellfun(@(x)(x.cueHalf.stats.bootPerf.p_percCorrect),lsrPerf_ch1);
[~,stats.cueHalf1.percCorrect.alpha_vsDelTow]    = FDR(stats.cueHalf1.percCorrect.pval_vsDelTow');
stats.cueHalf1.percCorrect.diffMean_CtrlvsDelTow = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch1);
stats.cueHalf1.percCorrect.diffStd_CtrlvsDelTow  = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.diff_percCorrect_std),lsrPerf_ch1);
stats.cueHalf1.percCorrect.pval_CtrlvsDelTow     = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.p_percCorrect),lsrPerf_ch1);
[~,stats.cueHalf1.percCorrect.alpha_CtrlvsDelTow]= FDR(stats.cueHalf1.percCorrect.pval_CtrlvsDelTow');

stats.cueHalf2.percCorrect.diffMean              = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch2);
stats.cueHalf2.percCorrect.diffStd               = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_std),lsrPerf_ch2);
stats.cueHalf2.percCorrect.pval                  = cellfun(@(x)(x.stats.bootPerf.p_percCorrect),lsrPerf_ch2);
[~,stats.cueHalf2.percCorrect.alpha]             = FDR(stats.cueHalf2.percCorrect.pval');
stats.cueHalf2.percCorrect.diffMean_vsDelTow     = cellfun(@(x)(x.cueHalf.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch2);
stats.cueHalf2.percCorrect.diffStd_vsDelTow      = cellfun(@(x)(x.cueHalf.stats.bootPerf.diff_percCorrect_std),lsrPerf_ch2);
stats.cueHalf2.percCorrect.pval_vsDelTow         = cellfun(@(x)(x.cueHalf.stats.bootPerf.p_percCorrect),lsrPerf_ch2);
[~,stats.cueHalf2.percCorrect.alpha_vsDelTow]    = FDR(stats.cueHalf2.percCorrect.pval_vsDelTow');
stats.cueHalf2.percCorrect.diffMean_CtrlvsDelTow = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch2);
stats.cueHalf2.percCorrect.diffStd_CtrlvsDelTow  = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch2);
stats.cueHalf2.percCorrect.pval_CtrlvsDelTow     = cellfun(@(x)(x.cueHalf_vsCtrl.stats.bootPerf.diff_percCorrect_mean),lsrPerf_ch2);
[~,stats.cueHalf2.percCorrect.alpha_CtrlvsDelTow]= FDR(stats.cueHalf2.percCorrect.pval_CtrlvsDelTow');

stats.delay.percCorrect.diffMean                 = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_mean),lsrPerf_del);
stats.delay.percCorrect.diffStd                  = cellfun(@(x)(x.stats.bootPerf.diff_percCorrect_std),lsrPerf_del);
stats.delay.percCorrect.pval                     = cellfun(@(x)(x.stats.bootPerf.p_percCorrect),lsrPerf_del);
[~,stats.delay.percCorrect.alpha]                = FDR(stats.delay.percCorrect.pval');

% logistic regression
logReg = {'decayIndex';'avgWeight'};
for iL = 1:numel(logReg)
  stats.cueHalf1.(['logReg_' logReg{iL}]).diffMean              = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_mean',logReg{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['logReg_' logReg{iL}]).diffStd               = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_std',logReg{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['logReg_' logReg{iL}]).pval                  = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_p',logReg{iL}))),lsrPerf_ch1);
  [~,stats.cueHalf1.(['logReg_' logReg{iL}]).alpha]             = FDR(stats.cueHalf1.(['logReg_' logReg{iL}]).pval');
   
  stats.cueHalf2.(['logReg_' logReg{iL}]).diffMean              = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_mean',logReg{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['logReg_' logReg{iL}]).diffStd               = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_std',logReg{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['logReg_' logReg{iL}]).pval                  = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_p',logReg{iL}))),lsrPerf_ch2);
  [~,stats.cueHalf2.(['logReg_' logReg{iL}]).alpha]             = FDR(stats.cueHalf2.(['logReg_' logReg{iL}]).pval');
  
  stats.delay.(['logReg_' logReg{iL}]).diffMean                 = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_mean',logReg{iL}))),lsrPerf_del);
  stats.delay.(['logReg_' logReg{iL}]).diffStd                  = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_diff_std',logReg{iL}))),lsrPerf_del);
  stats.delay.(['logReg_' logReg{iL}]).pval                     = cellfun(@(x)(eval(sprintf('x.lsr.logisticReg.%s_p',logReg{iL}))),lsrPerf_del);
  [~,stats.delay.(['logReg_' logReg{iL}]).alpha]                = FDR(stats.delay.(['logReg_' logReg{iL}]).pval');
end

% view angle
viewDecode = {'overall','cueHalf1','cueHalf2','mem'};
for iL = 1:numel(viewDecode)
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffMean              = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffStd               = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval                  = cellfun(@(x)(eval(sprintf('x.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch1);
  [~,stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).alpha]             = FDR(stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval');
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffMean_vsDelTow     = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffStd_vsDelTow      = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval_vsDelTow         = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch1);
  [~,stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).alpha_vsDelTow]    = FDR(stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval_vsDelTow');
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffMean_CtrlvsDelTow = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).diffStd_CtrlvsDelTow  = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch1);
  stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval_CtrlvsDelTow     = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch1);
  [~,stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).alpha_CtrlvsDelTow]= FDR(stats.cueHalf1.(['viewAngDecode_' viewDecode{iL}]).pval_CtrlvsDelTow');
  
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffMean              = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffStd               = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval                  = cellfun(@(x)(eval(sprintf('x.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch2);
  [~,stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).alpha]             = FDR(stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval');
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffMean_vsDelTow     = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffStd_vsDelTow      = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval_vsDelTow         = cellfun(@(x)(eval(sprintf('x.cueHalf.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch2);
  [~,stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).alpha_vsDelTow]    = FDR(stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval_vsDelTow');
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffMean_CtrlvsDelTow = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).diffStd_CtrlvsDelTow  = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_ch2);
  stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval_CtrlvsDelTow     = cellfun(@(x)(eval(sprintf('x.cueHalf_vsCtrl.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_ch2);
  [~,stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).alpha_CtrlvsDelTow]= FDR(stats.cueHalf2.(['viewAngDecode_' viewDecode{iL}]).pval_CtrlvsDelTow');
  
  stats.delay.(['viewAngDecode_' viewDecode{iL}]).diffMean                 = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_%s',viewDecode{iL}))),lsrPerf_del);
  stats.delay.(['viewAngDecode_' viewDecode{iL}]).diffStd                  = cellfun(@(x)(eval(sprintf('x.viewAngle.lsr.decodeAccDiff_sd_%s',viewDecode{iL}))),lsrPerf_del);
  stats.delay.(['viewAngDecode_' viewDecode{iL}]).pval                     = cellfun(@(x)(eval(sprintf('x.viewAngle.p_decodeAcc_%s',viewDecode{iL}))),lsrPerf_del);
  [~,stats.delay.(['viewAngDecode_' viewDecode{iL}]).alpha]                = FDR(stats.delay.(['viewAngDecode_' viewDecode{iL}]).pval');
end

end
