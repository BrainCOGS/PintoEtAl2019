function Pinto2019_figS2_ephys_varPower(ephysPath,varPowerPath,analysisFilePath,summaryFile)

% Pinto2019_figS2_ephys_varPower(ephysPath,varPowerPath,summaryFile)
% plots ephys and variable laser power behavioral data
% ephysPath is path for ephys data analysis files to be loaded
% varPowerPath is path for variable laser power behavioral data to be loaded
% analysisFilePath is path for laser data analysis files to be loaded (whole trial vgat, 2 or 6 mW)
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate ephys analysis files:
%    analyzeOptoRec_batch.m (repo: ephys) (will save summary_awake.mat)
%    requires sorted spikes (done here with klusta), refer to function for
%    more details
% 2) To generate variable lsr power behavioral analysis:
%    analyzeVarPower.m (will save the V1varPower.mat files)(repo: behav)
% 3) To generate inactivation experiment logs:
%    analyzeConcatLsrLog_batch(exptType). expType is a string that sets the
%    type of data (e.g. whole trial, subtrial etc). Relevant one for this
%    figure is 'wholeTrial'
%    (will save the fullGrid_*.mat files)(repo: behaviorAnalysis)
% 4) To statistically compare between inactivation experiments:
%    analyzeConcatLsrLog_xExptStats.m (will save
%    fullGrid_wholeTrial_xExptComp.mat)(repo: behaviorAnalysis)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.excitEgPath         = [ephysPath 'vg9awake_1/analysis.mat']; % rec path
cfg.excitEgAddress      = [3 4 1 1]; % consecutive indices in nested data structure
cfg.inhibitEgPath       = [ephysPath 'vg9awake_2/analysis.mat']; % rec path
cfg.inhibitEgAddress    = [4 4 1 3]; % consecutive indices in nested data structure
cfg.lsrDur              = 1.5; % sec
cfg.lsrFreq             = 40; % hz
cfg.rampDownDur         = 0.1; % sec
cfg.dists               = [0 .25 .5 1 2];
cfg.powers              = [.5 2 4 6 8];
cfg.ntrials             = 10;
cfg.sampRate            = 20000;
cfg.binSize             = 0.1;
cfg.baselineDur         = cfg.lsrDur;
cfg.reboundDur          = 0.5;
cfg.dutyCycle           = 0.8;
cfg.spkWdthCutoff       = 0.5;
cfg.depthBins           = 400:100:800;
cfg.reboundBins         = -100:20:100;
cfg.ipsiCl              = FormatDefaults.mediumGreen;
cfg.contraCl            = FormatDefaults.darkPurple;
cfg.V1gridID            = 26;
cfg.towersCl            = widefieldParams.darkgray; 
cfg.ctrlCl              = widefieldParams.darkgreen; 

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath = [getRepositoryPath(codeFile) '/stockFigs'];

%% load ephys opto data, compile summary stats
ephysSumm  = load([ephysPath 'summary_awake.mat']);
ephysStats = summaryStats_ephys(ephysSumm, cfg);
ephysSumm  = rmfield(ephysSumm,'data_all');

%% load behav var power data, compile summary stats
load(varPowerPath);
varPowerBehavStats = varPower.stats;

%% load accumulation maze data 2mw vs 6mw
thisdir = pwd; 
cd(analysisFilePath)
load fullGrid_lastMaze_whole_vgat_perfTh60_1mW_nBinLogReg4 lsrPerf
lsrPerf_1mw  = lsrPerf;
load fullGrid_lastMaze_whole_vgat_perfTh60_2mW_nBinLogReg4 lsrPerf
lsrPerf_2mw  = lsrPerf;
load fullGrid_lastMaze_whole_vgat_perfTh60_6mW_nBinLogReg4 lsrPerf
lsrPerf_6mw  = lsrPerf;
load fullGrid_wholeTrial_xExptComp pvals
xExptStats   = pvals;
cd(thisdir)

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'ephysSumm','ephysStats','varPower','lsrPerf_1mw','lsrPerf_2mw','lsrPerf_6mw','-v7.3')
  else
    save(summaryFile,'ephysSumm','ephysStats','varPower','lsrPerf_1mw','lsrPerf_2mw','lsrPerf_6mw','-append')
  end
end

%% Configure figure and panels
layout       = [ 1  1  2  2  4  6  8  ...
               ; 1  1  3  3  5  7  8  ...
               ; 9  10 11 12 13 14 15 ...
               ; 16 17 18 19 20 21 22 ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% Plot expt schematics
iPanel       = pnlExptScheme(fig, iPanel, cfg);

%% eg PSTHs + waveform
iPanel       = pnlPSTH(fig, iPanel, cfg);

%% effect vs power and distance (excit and inhibit)
iPanel       = pnlEffectVSpower(fig, iPanel, ephysStats, cfg);

%% effect vs depth
iPanel       = pnlDepth(fig, iPanel, ephysStats, cfg);

%% rebound histogram
iPanel       = pnlRebound(fig, iPanel, ephysStats, ephysSumm, cfg);

%% t4 schematic
iPanel       = pnlVisGuideSchem(fig, iPanel+1, cfg);

%% unilateral V1 schematic
iPanel       = pnlV1Schem(fig, iPanel, cfg);

%% perc correct + perc finish (ipsi vs contra)
iPanel       = pnlVarPower(fig, iPanel, varPower, cfg, 'vgat');

%% perc correct + perc finish (ipsi vs contra) WT
iPanel       = pnlVarPower(fig, iPanel, varPower, cfg, 'ctrl');

%% effect at different powers, towers task
iPanel       = pnlTaskComp_V1powers(fig, iPanel, lsrPerf_6mw, lsrPerf_2mw, lsrPerf_1mw, varPower, cfg);

%% main maze schematic
iPanel       = iPanel+2;
iPanel       = pnlMainMzSchem(fig, iPanel, cfg);

%% effect at different powers, towers task
[iPanel,stats] = pnlPowerComp(fig, iPanel, lsrPerf_6mw, lsrPerf_2mw, lsrPerf_1mw);

%% Save figure to disk
versionInfo.ephysStats         = ephysStats;
versionInfo.varPowerBehavStats = varPowerBehavStats;
versionInfo.powerCompStats     = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% summary stats
function stats = summaryStats_ephys(ephysSumm, cfg)

% n
stats.nPenetrations  = numel(ephysSumm.awake_recls);
stats.mice           = unique(cellfun(@(x)(x(1:3)),ephysSumm.awake_recls,'UniformOutput',false));
stats.nmice          = numel(stats.mice);
stats.nUnits         = numel(ephysSumm.isSUA);
stats.nSUA           = sum(ephysSumm.isSUA);
stats.nMUA           = sum(~ephysSumm.isSUA);
isFS                 = ephysSumm.spkWidth < cfg.spkWdthCutoff;
stats.isFS           = isFS;
stats.nFS            = sum(stats.isFS);
stats.nRS            = sum(~stats.isFS);

% mean, std effects
stats.summMatCols    = 'distance';
stats.summMatRows    = 'power';
stats.summMatFS_mean = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatFS_std  = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatRS_mean = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatRS_std  = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatFS_rebound_mean = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatFS_rebound_std  = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatRS_rebound_mean = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatRS_rebound_std  = zeros(numel(cfg.powers),numel(cfg.dists));
stats.summMatRS_depth_mean   = zeros(numel(cfg.powers),numel(cfg.depthBins)-1);
stats.summMatRS_depth_std    = zeros(numel(cfg.powers),numel(cfg.depthBins)-1);
stats.summMatRS_depth_mean   = zeros(numel(cfg.powers),numel(cfg.depthBins)-1);
stats.summMatRS_depth_std    = zeros(numel(cfg.powers),numel(cfg.depthBins)-1);

for iP = 1:numel(cfg.powers)
  for iD = 1:numel(cfg.dists)
    stats.summMatFS_mean(iP,iD) = 100.*nanmean(ephysSumm.supprMat(isFS,iD,iP));
    stats.summMatFS_std(iP,iD)  = 100.*nanstd(ephysSumm.supprMat(isFS,iD,iP))./sqrt(stats.nFS-1);
    stats.summMatRS_mean(iP,iD) = 100.*nanmean(ephysSumm.supprMat(~isFS,iD,iP));
    stats.summMatRS_std(iP,iD)  = 100.*nanstd(ephysSumm.supprMat(~isFS,iD,iP))./sqrt(stats.nRS-1);
    stats.summMatFS_rebound_mean(iP,iD) = 100.*nanmean(ephysSumm.reboundMat(isFS,iD,iP));
    stats.summMatFS_rebound_std(iP,iD)  = 100.*nanstd(ephysSumm.reboundMat(isFS,iD,iP))./sqrt(stats.nFS-1);
    stats.summMatRS_rebound_mean(iP,iD) = 100.*nanmean(ephysSumm.reboundMat(~isFS,iD,iP));
    stats.summMatRS_rebound_std(iP,iD)  = 100.*nanstd(ephysSumm.reboundMat(~isFS,iD,iP))./sqrt(stats.nRS-1);

    % by depth at dist = 0
    if iD > 1; continue; end
    for iB = 1:numel(cfg.depthBins)-1
      isDepth                           = ephysSumm.depth > cfg.depthBins(iB) & ephysSumm.depth <= cfg.depthBins(iB+1);
      stats.summMatFS_depth_mean(iP,iB) = 100.*nanmean(ephysSumm.supprMat(isFS & isDepth,iD,iP));
      stats.summMatFS_depth_sem(iP,iB)  = 100.*nanstd(ephysSumm.supprMat(isFS & isDepth,iD,iP))./sqrt(stats.nFS-1);
      stats.summMatRS_depth_mean(iP,iB) = 100.*nanmean(ephysSumm.supprMat(~isFS & isDepth,iD,iP));
      stats.summMatRS_depth_sem(iP,iB)  = 100.*nanstd(ephysSumm.supprMat(~isFS & isDepth,iD,iP))./sqrt(stats.nRS-1);
    end
    
  end
end

end

%% ------------------------------------------------------------------------
%% expt schematics
function iPanel = pnlExptScheme(fig,iPanel,cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/ephysSchematics.png']);
imshow(im); axis off; axis image

end

%% ------------------------------------------------------------------------
%% eg PSTHs
function iPanel  = pnlPSTH(fig, iPanel, cfg)

%% excitatory eg
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

%% laser waveform
rgRate = 400;
lsrVec = [];
szl    = ((1000*cfg.lsrDur)/cfg.lsrFreq)*(rgRate/(1000*cfg.lsrDur));
lsrVec = [lsrVec; repmat([ones(floor(szl*cfg.dutyCycle),1); ...
          zeros(ceil(szl*(1-cfg.dutyCycle)),1)],cfg.lsrFreq,1)];
lsrVec = repmat(lsrVec,[ceil(cfg.lsrDur) 1]);
lsrVec = lsrVec(1:round(cfg.lsrDur*rgRate));
rampD  = [ones(rgRate*(cfg.lsrDur-cfg.rampDownDur),1); linspace(1,0,rgRate*cfg.rampDownDur)'];
lsrVec = [zeros(rgRate*cfg.baselineDur,1); lsrVec.*rampD; zeros(rgRate*cfg.baselineDur,1)];

%% plot psth with lsr and waveform inset 
load(cfg.excitEgPath,'data')
idx      = cfg.excitEgAddress;
thiscell = data.loc(idx(1)).pow(idx(2)).dur(idx(3)).clust(idx(4));
wvf      = data.waveforms.avg(idx(4),:);
psth_t   = thiscell.psth_t;
tidx     = find(psth_t >= -cfg.baselineDur,1,'first'):find(psth_t <= cfg.baselineDur*2,1,'last');
psth_t   = psth_t(tidx);

%%
bar(psth_t,mean(thiscell.psth(:,tidx)),...
    'facecolor',FormatDefaults.mediumGray,'edgecolor',FormatDefaults.mediumGray,'barwidth',1);
yl  = get(gca,'ylim');
lsr = yl(2) + diff(yl)*.2*lsrVec;
plot(linspace(psth_t(1),psth_t(end),numel(lsr)),lsr,'-','color',analysisParams.lsrCl);
text(-cfg.baselineDur+.1,min(lsr)+.5,'Laser','color',analysisParams.lsrCl,'fontsize',9);
wvx = linspace(cfg.baselineDur*2-.5,cfg.baselineDur*2,numel(wvf));
wva = max(wvf)-min(wvf);
wvr = diff(yl)*.18;
wvt = wvf * (wvr/wva) + yl(2) + diff(yl)*.2;
plot(wvx,wvt,'-','color',FormatDefaults.darkGray,'linewidth',.75)
plot([wvx(1)-.1 wvx(cfg.sampRate/1000)-.1],[yl(2)+diff(yl)*.03 yl(2)+diff(yl)*.03],'k-','linewidth',.25)
plot([wvx(1)-.1 wvx(1)-.1],[yl(2)+diff(yl)*.03 yl(2)+diff(yl)*.03+wvr/2],'k-','linewidth',.25)
text(wvx(1)-.09,yl(2)+diff(yl)*.013,'1 ms','fontsize',8,'horizontalAlignment','left')
th  = text(wvx(1)-.14,yl(2)+diff(yl)*.035,'50 \muV','fontsize',8,'horizontalAlignment','left');
set(th,'rotation',90)

xlabel('Time from laser onset (s)')
ylabel('FR (Hz)')
set(gca,'ytick',0:5:5*round(yl(2)/5));
xlim([-cfg.baselineDur cfg.baselineDur*2])

%% inhibitory eg
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

%% plot psth with lsr and waveform inset 
load(cfg.inhibitEgPath,'data')
idx      = cfg.inhibitEgAddress;
thiscell = data.loc(idx(1)).pow(idx(2)).dur(idx(3)).clust(idx(4));
wvf      = data.waveforms.avg(idx(4),:);
psth_t   = thiscell.psth_t;
tidx     = find(psth_t >= -cfg.baselineDur,1,'first'):find(psth_t <= cfg.baselineDur*2,1,'last');
psth_t   = psth_t(tidx);

%%
bar(psth_t,mean(thiscell.psth(:,tidx)),...
    'facecolor',FormatDefaults.mediumGray,'edgecolor',FormatDefaults.mediumGray,'barwidth',1);
yl  = get(gca,'ylim');
lsr = yl(2) + diff(yl)*.2*lsrVec;
plot(linspace(psth_t(1),psth_t(end),numel(lsr)),lsr,'-','color',analysisParams.lsrCl);
text(-cfg.baselineDur+.1,min(lsr)+.5,'Laser','color',analysisParams.lsrCl,'fontsize',9);
wvx = linspace(cfg.baselineDur*2-.5,cfg.baselineDur*2,numel(wvf));
wva = max(wvf)-min(wvf);
wvr = diff(yl)*.18;
wvt = wvf * (wvr/wva) + yl(2) + diff(yl)*.2;
plot(wvx,wvt,'-','color',FormatDefaults.darkGray,'linewidth',.75)
plot([wvx(1)-.1 wvx(cfg.sampRate/1000)-.1],[yl(2)+diff(yl)*.03 yl(2)+diff(yl)*.03],'k-','linewidth',.25)
plot([wvx(1)-.1 wvx(1)-.1],[yl(2)+diff(yl)*.03 yl(2)+diff(yl)*.03+wvr/2],'k-','linewidth',.25)
text(wvx(1)-.09,yl(2)+diff(yl)*.013,'1 ms','fontsize',8,'horizontalAlignment','left')
th  = text(wvx(1)-.14,yl(2)+diff(yl)*.035,'50 \muV','fontsize',8,'horizontalAlignment','left');
set(th,'rotation',90)

xlabel('Time from laser onset (s)')
ylabel('FR (Hz)')
set(axs,'ytick',0:50:50*round(yl(2)/50));
xlim([-cfg.baselineDur cfg.baselineDur*2])

end

%% ------------------------------------------------------------------------
%% effect vs power vs lateral distance
function iPanel= pnlEffectVSpower(fig, iPanel, ephysStats, cfg)

cmap = fliplr(red(numel(cfg.powers)*2+1));
cmap = cmap(2:2:end,:);

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

for iP = 1:numel(cfg.powers)
  errorbar(cfg.dists, ephysStats.summMatRS_mean(iP,:), ephysStats.summMatRS_std(iP,:), ...
           '-', 'color', cmap(iP,:))
end
set(axs, 'xtick', cfg.dists)

xlabel('Lateral dist. from beam (mm)')
ylabel(sprintf('Suppression (%%)\n(Putative excit.)'))
xlim([-.05 2.05]); ylim([0 100])

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

for iP = 1:numel(cfg.powers)
  errorbar(cfg.dists, -ephysStats.summMatFS_mean(iP,:), ephysStats.summMatRS_std(iP,:), ...
           '-', 'color', cmap(iP,:))
end
legend(cellfun(@(x)([num2str(x) ' mW']),num2cell(cfg.powers),'UniformOutput',false),'location','best','fontsize',8)
set(axs, 'xtick', cfg.dists)

xlabel('Lateral dist. from beam (mm)')
ylabel(sprintf('Excitation (%%)\n(Putative FS.)'))
xlim([-.05 2.05]); ylim([0 4500])

end

%% ------------------------------------------------------------------------
%% effect vs power vs depth
function iPanel= pnlDepth(fig, iPanel, ephysStats, cfg)

cmap = fliplr(red(numel(cfg.powers)*2+1));
cmap = cmap(2:2:end,:);

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

for iP = 1:numel(cfg.powers)
  errorbar(toBinCenters(cfg.depthBins), ephysStats.summMatRS_depth_mean(iP,:), ...
           ephysStats.summMatRS_depth_std(iP,:), '-', 'color', cmap(iP,:))
end

set(axs, 'xtick', 200:200:800, 'ytick', 0:20:100)
xlim([400 800]); ylim([0 100])

xlabel('Depth from surface (mm)')
ylabel(sprintf('Suppression (%%)\n(Putative excit.)'))

set(axs,'view',[90 90])

end

%% ------------------------------------------------------------------------
%% rebound
function iPanel = pnlRebound(fig, iPanel, ephysStats, ephysSumm, cfg)

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

dist = cfg.dists == 1;
pow  = cfg.powers == 6;
reb  = 100.*ephysSumm.reboundMat(~ephysStats.isFS,dist,pow);

plot([0 0],[0 20],'--','color',FormatDefaults.lightGray)
bins                    = cfg.reboundBins;
xaxis                   = toBinCenters(bins);
counts                  = histcounts(reb,bins);
bar ( axs, xaxis, counts, FormatDefaults.barWidth                 ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , FormatDefaults.mediumGray                 ...
                              );

% plot arrowheads to indicate mean
plot(axs, nanmedian(reb), .6,  'v',                                             ...
                             'markerfacecolor',     'w',                                    ...
                             'markeredgecolor',     FormatDefaults.darkGray,                ...
                             'linewidth'      ,     FormatDefaults.linewidthThick)                            

set(axs, 'xtick', -100:25:100, 'xticklabel',{'-100','','-50','','0','','50','','100'})
xlim([-100 100]); ylim([0 20])

xlabel('Rebound (%)')
ylabel('Num. RS units')

end

%% ------------------------------------------------------------------------
%% visual guide schematics
function iPanel = pnlVisGuideSchem(fig, iPanel, cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematicVisualGuideNoDistr.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-20,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)
end

%% ------------------------------------------------------------------------
%% expt schematics V1
function iPanel = pnlV1Schem(fig, iPanel, cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/allenBrainOverlay.tiff']);
im     = repmat(im(:,:,1),[1 1 3]);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
bregma   = [285 274];
PxlPerMM = 41.4;
data     = load([analysisParams.gridpath 'V1unilateralSingleGrid.mat'],'grid'); 
grid     = data.grid;

for iSide = 1:size(grid,1)
    x = grid(iSide,1)*PxlPerMM + bregma(2);
    y = grid(iSide,2)*PxlPerMM + bregma(1);
    if iSide == 1
      plot(x,y,'o','color',analysisParams.lsrShade,'markerfacecolor',analysisParams.lsrCl,'markersize',4)
    else
      plot(x,y,'x','color',analysisParams.lsrShade,'linewidth',.5,'markersize',4)
    end
end
end

%% ------------------------------------------------------------------------
%% behav var power, % correct and excess travel
function iPanel = pnlVarPower(fig, iPanel, varPower, cfg, plotWhat)

% % correct
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

powers          = varPower.powers;
pc_ctrl         = varPower.(plotWhat).percCorrect_ctrl(1);
pc_ctrl_seml    = varPower.(plotWhat).percCorrect_ctrl_binoErr(1,1)*100;
pc_ctrl_semu    = varPower.(plotWhat).percCorrect_ctrl_binoErr(1,2)*100;

pc_ipsi         = varPower.(plotWhat).percCorrect_ipsi';
pc_ipsi_seml    = -varPower.(plotWhat).percCorrect_ipsi_binoErr(:,1).*100 + pc_ipsi;
pc_ipsi_semu    = varPower.(plotWhat).percCorrect_ipsi_binoErr(:,2).*100 - pc_ipsi;

pc_contra       = varPower.(plotWhat).percCorrect_contra';
pc_contra_seml  = -varPower.(plotWhat).percCorrect_contra_binoErr(:,1).*100 + pc_contra;
pc_contra_semu  = varPower.(plotWhat).percCorrect_contra_binoErr(:,2).*100 - pc_contra;

h1 = plot([powers(1) powers(end)],[pc_ctrl pc_ctrl],'-','color',analysisParams.ctrlCl);
plot([powers(1) powers(end)],[pc_ctrl_seml pc_ctrl_seml],'--','color',analysisParams.ctrlCl,'linewidth',.25);
plot([powers(1) powers(end)],[pc_ctrl_semu pc_ctrl_semu],'--','color',analysisParams.ctrlCl,'linewidth',.25);

h2 = errorbar(powers,pc_ipsi,pc_ipsi_seml,pc_ipsi_semu,'.-','markersize',10,'color',cfg.ipsiCl);
h3 = errorbar(powers,pc_contra,pc_contra_seml,pc_contra_semu,'.-','markersize',10,'color',cfg.contraCl);

ylabel('Perf. (% correct)'); 
xlabel('Laser power (mW)');
set(axs,'xscale','log','xtick',powers); 
xlim([.21 15]); ylim([0 100]);
switch plotWhat
  case 'vgat'
    title('ChR2+ (vgat)','fontsize',13,'color','k')
  case 'ctrl'
    title('ChR2- (wt)','fontsize',13,'color','k')
    legend([h1 h2 h3],{'laser off','ipsi','contra'},'location','best','fontsize',8)
end

end

%% ------------------------------------------------------------------------
%% visual guide schematics
function iPanel = pnlMainMzSchem(fig, iPanel, cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
plot([-10 -10],[1 nX],'-','color',analysisParams.lsrCl,'linewidth',3)
t = text(-20,(nX-1)/2,'Laser on','color',analysisParams.lsrCl,...
        'horizontalAlignment','center','fontsize',FormatDefaults.legendFontSize);
set(t,'rotation',90)
end

%% ------------------------------------------------------------------------
%% expt schematics V1
function iPanel = pnlGridSchem(fig, iPanel, cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/allenBrainOverlay.tiff']);
im     = repmat(im(:,:,1),[1 1 3]);

hold(axs, 'on')
image(flipud(im)); axis off; axis image
bregma   = [285 274];
PxlPerMM = 41.4;
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
%% delta performance, by power, with stats
function [iPanel,stats] = pnlPowerComp(fig, iPanel, lsrPerf_6mw, lsrPerf_2mw, lsrPerf_1mw)

%% effect sizes
[pvals6,alphac6]               = retrievePvals(lsrPerf_6mw,'percCorrect','boot');
stats.mainMaze6mw.perf.alpha   = alphac6;
stats.mainMaze6mw.perf.pvals   = pvals6.p(1:2:end);
stats.mainMaze6mw.perf.effSize = pvals6.effectSize(1:2:end);
isSig6mw                       = pvals6.p(1:2:end) <= alphac6;
 
[pvals2,alphac2]               = retrievePvals(lsrPerf_2mw,'percCorrect','boot');
stats.mainMaze2mw.perf.alpha   = alphac2;
stats.mainMaze2mw.perf.pvals   = pvals2.p(1:2:end);
stats.mainMaze2mw.perf.effSize = pvals2.effectSize(1:2:end);
isSig2mw                       = pvals2.p(1:2:end) <= alphac2;
 
[pvals1,alphac1]               = retrievePvals(lsrPerf_1mw,'percCorrect','boot');
stats.mainMaze1mw.perf.alpha   = alphac1;
stats.mainMaze1mw.perf.pvals   = pvals1.p(1:2:end);
stats.mainMaze1mw.perf.effSize = pvals1.effectSize(1:2:end);
isSig1mw                       = pvals1.p(1:2:end) <= alphac1;
 
effSize = 100.*[stats.mainMaze1mw.perf.effSize stats.mainMaze2mw.perf.effSize stats.mainMaze6mw.perf.effSize];

%% ANOVA
stats.mainMaze.powers.lbls            = {'1 mW','2 mW','6mW'};
stats.mainMaze.powers.effectSize_mean = [mean(stats.mainMaze1mw.perf.effSize) mean(stats.mainMaze2mw.perf.effSize) mean(stats.mainMaze6mw.perf.effSize)];
stats.mainMaze.powers.effectSize_sem  = [std(stats.mainMaze1mw.perf.effSize) std(stats.mainMaze2mw.perf.effSize) std(stats.mainMaze6mw.perf.effSize)]./sqrt(numel(isSig2mw)-1);
stats.mainMaze.powers.nSig            = [sum(isSig1mw) sum(isSig2mw) sum(isSig6mw)];

powersVec = [ones(size(isSig1mw)) ones(size(isSig1mw))+1 ones(size(isSig1mw))+2];
roiVec    = [(1:numel(isSig1mw))' (1:numel(isSig1mw))'   (1:numel(isSig1mw))'  ];
[stats.mainMaze.powers.ANOVA_pvals,stats.mainMaze.powers.ANOVA_table,stats.mainMaze.powers.ANOVA_stats] ...
          = anovan(effSize(:),{powersVec(:),roiVec(:)},'varnames',{'power','location'},'display','off');


%% plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on') 
plotBrainPVals(pvals1.coord(:,1),pvals1.coord(:,2),pvals1.p,pvals1.effectSize,alphac1,1,'es_size',axs,1);
title('\Delta Perf. (%) 1 mW','fontsize',13,'fontweight','bold')
 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on') 
plotBrainPVals(pvals2.coord(:,1),pvals2.coord(:,2),pvals2.p,pvals2.effectSize,alphac2,1,'es_size',axs,0);
title('\Delta Perf. (%) 2 mW','fontsize',13,'fontweight','bold')
 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on') 
plotBrainPVals(pvals6.coord(:,1),pvals6.coord(:,2),pvals6.p,pvals6.effectSize,alphac6,1,'es_size',axs,0);
title('\Delta Perf. (%) 6 mW','fontsize',13,'fontweight','bold')
 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on') 
plot(repmat([1 2 6],[size(effSize,1) 1])', effSize', '-', 'color', FormatDefaults.lightGray)
plot(1, effSize(~isSig1mw,1), 'o', 'color', FormatDefaults.lightGray, 'markerfacecolor', FormatDefaults.lightGray,'markersize',3);
plot(2, effSize(~isSig2mw,2), 'o', 'color', FormatDefaults.lightGray, 'markerfacecolor', FormatDefaults.lightGray,'markersize',3);
plot(1, effSize(isSig1mw,1), 'o', 'color', FormatDefaults.darkGray, 'markerfacecolor', FormatDefaults.darkGray,'markersize',3);
plot(2, effSize(isSig2mw,2), 'o', 'color', FormatDefaults.darkGray, 'markerfacecolor', FormatDefaults.darkGray,'markersize',3);
if sum(~isSig6mw)>0
  plot(6, effSize(~isSig6mw,3), 'o', 'color', FormatDefaults.lightGray, 'markerfacecolor', FormatDefaults.lightGray,'markersize',3);
end
plot(6, effSize(isSig6mw,3), 'o', 'color', FormatDefaults.darkGray, 'markerfacecolor', FormatDefaults.darkGray,'markersize',3);
 
xlim([.9 6.6]); ylim([-30 0])
text(3,0,sprintf('p(power) = %1.2g\np(region) = %1.2g',stats.mainMaze.powers.ANOVA_pvals(1),stats.mainMaze.powers.ANOVA_pvals(2)),...
     'horizontalAlignment','left','color','k','fontsize',9);
set(axs,'xtick', [1 2 6])
ylabel('\Delta perf. (%) by location'); xlabel('Laser power (mW)')
set(axs,'xscale','log')

end


%% compare effect of varying power in V1 across different tasks
function iPanel = pnlTaskComp_V1powers(fig, iPanel, lsrPerf_6mw, lsrPerf_2mw, lsrPerf_1mw, varPower, cfg)

%% get effect sizes
% normalize such that control perf is 100% and chance is 0%
v1idx                = cfg.V1gridID;
accumul_perfDrop     = [];
accumul_binoErr      = [];
accumul_perfDrop(1)  = 0;
accumul_binoErr(1,:) = 0;

accumul_perfDrop(end+1)  = abs(lsrPerf_1mw{v1idx}.lsr.percCorrect-lsrPerf_1mw{v1idx}.ctrl.percCorrect) / (lsrPerf_1mw{v1idx}.ctrl.percCorrect - 50);
accumul_binoErr(end+1,:) = (lsrPerf_1mw{v1idx}.stats.bootPerf.diff_percCorrect_std) ./ (lsrPerf_1mw{v1idx}.ctrl.percCorrect - 50);

accumul_perfDrop(end+1)  = abs(lsrPerf_2mw{v1idx}.lsr.percCorrect-lsrPerf_2mw{v1idx}.ctrl.percCorrect) / (lsrPerf_2mw{v1idx}.ctrl.percCorrect - 50);
accumul_binoErr(end+1,:) = (lsrPerf_2mw{v1idx}.stats.bootPerf.diff_percCorrect_std) ./ (lsrPerf_2mw{v1idx}.ctrl.percCorrect - 50);

accumul_perfDrop(end+1)  = abs(lsrPerf_6mw{v1idx}.lsr.percCorrect-lsrPerf_6mw{v1idx}.ctrl.percCorrect) / (lsrPerf_6mw{v1idx}.ctrl.percCorrect - 50);
accumul_binoErr(end+1,:) = (lsrPerf_6mw{v1idx}.stats.bootPerf.diff_percCorrect_std) ./ (lsrPerf_6mw{v1idx}.ctrl.percCorrect - 50);
accumul_x                = [0 1 2 6];
 

visguide_perfdrop = abs( varPower.vgat.percCorrect_both - varPower.vgat.percCorrect_ctrl) ./ (varPower.vgat.percCorrect_ctrl - 50);
visguide_binoErr  = (varPower.stats.vgat.sdDiff_percCorrect_both) ./ (varPower.vgat.percCorrect_ctrl./100 - 0.5);
visguide_x        = varPower.powers;
keepidx           = visguide_x == 1 | visguide_x == 2;
avgidx            = visguide_x == 4 | visguide_x == 8;
visguide_perfdrop = [0 visguide_perfdrop(keepidx) mean(visguide_perfdrop(avgidx))];
visguide_binoErr  = [0 visguide_binoErr(keepidx)  mean(visguide_binoErr(avgidx))];

%% plot
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on') 

plot([0 6],[0 0],'--','color',[.7 .7 .7])
plot([0 6],[100 100],'--','color',[.7 .7 .7])
h(1) = errorbar(accumul_x,100 - 100*accumul_perfDrop,100*accumul_binoErr','o-','color',cfg.towersCl,'markerfacecolor',cfg.towersCl,'markersize',3);
h(2) = errorbar(accumul_x,100 - 100*visguide_perfdrop,100*visguide_binoErr,'o-','color',cfg.ctrlCl,'markerfacecolor',cfg.ctrlCl,'markersize',3);
xlabel('Laser power (mW)'); 
ylabel('Performance'); 
text(3,95,'Laser off performance','color',[.7 .7 .7])
text(3,-5,'Chance performance','color',[.7 .7 .7])
set(gca,'ytick',0:25:100,'xtick',accumul_x); 
% ylim([-20 120])
xlim([-.5 6.5])
title('V1 inactivation')
legend(h,{'towers','vis guided'},'location','best')
end