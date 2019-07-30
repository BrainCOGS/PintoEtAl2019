function Pinto2019_figS2_wholeTrialOpto_WTctrl(analysisFilePath)

% Pinto2019_figS2_wholeTrialOpto_WTctrl(analysisFilePath)
% plots whole-trial inactivation for wt controls
% analysisFilePath is path for data analysis files to be loaded

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate inactivation experiment logs:
%    analyzeConcatLsrLog_batch(exptType). expType is a string that sets the
%    type of data (e.g. whole trial, subtrial etc). Relevant ones for this
%    figure are 'wholeTrial', 'WTctrls'
%    (will save the fullGrid_*.mat files)(repo: behaviorAnalysis)
% 2) To statistically compare between inactivation experiments:
%    analyzeConcatLsrLog_xExptStats.m (will save
%    fullGrid_wholeTrial_xExptComp.mat)(repo: behaviorAnalysis)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.perfTh              = 0.6;
cfg.nBins               = 4;
cfg.logRegBins          = linspace(10,200,cfg.nBins+1);
cfg.WTplotWhat          = {'percCorrect','bias_abs','viewDecode','speed','excessTravel'};
cfg.WTplotWhatLbl       = {'perf. (%)','|bias| (%)','unusualEvents','speed (%)','excess travel (%)'};

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

%% load opto data
thisdir = pwd; 
cd(analysisFilePath)
load fullGrid_lastMaze_whole_wt_perfTh60_nBinLogReg4 lsrPerf
lsrPerf_wt   = lsrPerf;
load fullGrid_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf
load fullGrid_wholeTrial_xExptComp pvals
xExptStats_fullGrid = pvals;
cd(thisdir)

%% SUPPL FIGURE: wt ctrl maps 
stats             = [];
layout            = [1 2 4 6 ...
                    ;1 3 5 7 ...
                    ];
fig              = PaneledFigure(layout, 'smaller');
stats             = figWT(fig, lsrPerf_wt, lsrPerf, cfg, stats, xExptStats_fullGrid);
versionInfo.stats = stats;

fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% wt ctrl maps
function stats = figWT(fig, lsrPerf, lsrPerfvgat, cfg, stats, xExptStats)

% maze schematic 
iPanel = 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-30,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)

% wt vs vgat
[pvals,alpha] = retrievePvals(lsrPerf,'percCorrect','boot');
effSizewt     = pvals.effectSize(1:2:end);
isSigwt       = pvals.p(1:2:end) <= alpha;
[pvals,alpha] = retrievePvals(lsrPerfvgat,'percCorrect','boot');
effSizevg     = pvals.effectSize(1:2:end);
isSigvg       = pvals.p(1:2:end) <= alpha;

effSize = 100.*[effSizewt effSizevg];

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

plot([.5 2.5], [0 0], '--', 'color', FormatDefaults.lightGray, 'linewidth', .25)
plot(repmat([1 2],[size(effSize,1) 1])', effSize', '-', 'color', FormatDefaults.lightGray)
h1 = plot(1, effSize(~isSigwt,1), 'o', 'color', FormatDefaults.lightGray, 'markerfacecolor', FormatDefaults.lightGray);
if sum(isSigwt)>0;  plot(1, effSize(isSigwt,1), 'o', 'color', FormatDefaults.darkGray, 'markerfacecolor', FormatDefaults.darkGray);    end
if sum(~isSigvg)>0; plot(2, effSize(~isSigvg,2), 'o', 'color', FormatDefaults.lightGray, 'markerfacecolor', FormatDefaults.lightGray); end
h2 = plot(2, effSize(isSigvg,2), 'o', 'color', FormatDefaults.darkGray, 'markerfacecolor', FormatDefaults.darkGray);
h3 = plot(2, effSize(xExptStats.isSig_vgatvsCtrl,2), 'ro', 'linewidth', .25);
legend([h1(1); h2(1); h3(1)], {'not signif.','signif.','signif. vs WT'}, 'location', 'best');
xlim([.5 2.5])
set(axs,'xtick', 1:2, 'xticklabel', {'ChR2- (WT)','ChR2+ (VGAT)'})
ylabel('\Delta perf. (%) by location')

% maps
for iP = 1:numel(cfg.WTplotWhat)
  iPanel = iPanel+1;
  axs    = fig.panel(iPanel);
  [pvals,alphac]  = retrievePvals(lsrPerf,cfg.WTplotWhat{iP},'boot');
  [~,cbar] = plotBrainPVals(pvals.coord(:,1),pvals.coord(:,2),pvals.p,pvals.effectSize,alphac,1,'es_size',axs,iP==1);
  title(['\Delta ' cfg.WTplotWhatLbl{iP}],'fontsize',12)
  
  stats.mainMazeWT.(cfg.WTplotWhat{iP}).alpha   = alphac;
  stats.mainMazeWT.(cfg.WTplotWhat{iP}).pvals   = pvals.p(1:2:end);
  stats.mainMazeWT.(cfg.WTplotWhat{iP}).effSize = pvals.effectSize(1:2:end);
  
  if iP == 1 
    pos             = cbar.Position;
    cbar.Units      = get(axs,'units');
    cbar.Position   = pos;
  end
end

end