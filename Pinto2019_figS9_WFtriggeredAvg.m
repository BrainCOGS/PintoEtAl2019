function Pinto2019_figS9_WFtriggeredAvg(analysisFilePath,summaryFile)

% Pinto2019_figS9_WFtriggeredAvg(analysisFilePath,summaryFile)
% plots event-triggered averages of df/f
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) see owf_WFdynamics.m for preprocessing
% 2) To generate triggered averages for each rec: (repo: widefieldImaging)
%    avgDff_eventTriggered.m (use default inputs)
% 3) To summarize activity averages: (repo: widefieldImaging)
%    summarizeAvgDffTrig.m (will save avgDffSummary.mat, loaded here)
%% ------------------------------------------------------------------------

%% analysis config
cfg.events          = {'tower','turn','preturn'};
cfg.egrecs          = {'ai3/20170126','ai7/20170814','ai10/20180409'};
cfg.egFrameIdx      = 4:2:20;
nF                  = numel(cfg.egFrameIdx);
cfg.egFrameIdx_pt   = 7:2:5+2*nF;
cfg.egTrialType     = 'correct';
cfg.egSide          = 'L';
cfg.colorMap        = parula;
cfg.compressPrctile = [1 99];

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
cd(analysisFilePath)
avgTrigDff            = load('avgDffSummary_triggered.mat');
avgTrigDff.avgDffSumm = rmfield(avgTrigDff.avgDffSumm,'mouse');

for iRec = 1:numel(cfg.egrecs)
  cd([analysisFilePath cfg.egrecs{iRec}])
  load avgDffEventTrig avgDff 
  load ROIfromRef ROI 
  
  avgTrigDff.mapEg.(cfg.events{iRec}).ROI = ROI;
  
  switch cfg.events{iRec}
    case 'preturn'
      avgTrigDff.mapEg.(cfg.events{iRec}).frames = avgDff.maze(end).([cfg.egTrialType cfg.egSide]).([cfg.events{iRec} cfg.egSide]).avg(:,:,cfg.egFrameIdx_pt);
      avgTrigDff.mapEg.(cfg.events{iRec}).taxis  = avgDff.timeAxis_preturb(cfg.egFrameIdx_pt);
    otherwise
      avgTrigDff.mapEg.(cfg.events{iRec}).frames = avgDff.maze(end).([cfg.egTrialType cfg.egSide]).([cfg.events{iRec} cfg.egSide]).avg(:,:,cfg.egFrameIdx);
      avgTrigDff.mapEg.(cfg.events{iRec}).taxis  = avgDff.timeAxis(cfg.egFrameIdx);
  end
end

cd(analysisFilePath)

%% Configure figure and panels
layout       = [ 1:nF        ...
               ; nF+1:nF*2   ...
               ; nF*2+1:nF*3 ...
               ; nF*3+1:nF*4 ...
               ];
fig          = PaneledFigure(layout, 'tiny');
iPanel       = 0;

%% frames (rec egs)
iPanel       = pnlFrames(fig,iPanel,avgTrigDff.mapEg.tower,cfg);
iPanel       = pnlFrames(fig,iPanel,avgTrigDff.mapEg.turn,cfg);
iPanel       = pnlFrames(fig,iPanel,avgTrigDff.mapEg.preturn,cfg);

%% averages
iPanel       = pnlROIavg(fig,iPanel,avgTrigDff.avgDffSumm,cfg);

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'avgTrigDff','-v7.3')
  else
    save(summaryFile,'avgTrigDff','-append')
  end
end

%% Save figure to disk
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% frames
function iPanel = pnlFrames(fig,iPanel,data,cfg)

cmap      = cfg.colorMap;
cmap(1,:) = [0 0 0];
clim      = [prctile(data.frames(:),cfg.compressPrctile(1)) prctile(data.frames(:),cfg.compressPrctile(2))];
barsize   = widefieldParams.pxlPerMM / widefieldParams.dsFactor;
nF        = size(data.frames,3);

for iF = 1:nF
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  hold(axs,'on')
  imagesc(data.frames(:,:,iF),clim);
  drawROI(data.ROI);
  colormap(cmap); axis image; axis off; axis ij
  text(5,140,sprintf('%1.1f s',data.taxis(iF)),'color','k')
  plot([5 5+barsize],[5 5],'w-','linewidth',2.5); % scale bar
end

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\DeltaF/F (z-score)';

end

%% ------------------------------------------------------------------------
%% ROI averages
function iPanel = pnlROIavg(fig,iPanel,avgDffSumm,cfg)

lblsAll   = avgDffSumm.ROIlbl;
lbls      = unique(cellfun(@(x)(x(1:end-2)),lblsAll,'UniformOutput',false),'stable');
[cl,lbls] = getDefaultROIcl(lbls);
nROI      = numel(lblsAll);
avgAT     = [];

for iEv = 1:numel(cfg.events)
  avgAT(:,:,1)   = nanmean(avgDffSumm.(['ROIavg_accumul_' (cfg.egTrialType)]).([cfg.events{iEv} 'R']),3);
  avgAT(:,:,2)   = nanmean(avgDffSumm.(['ROIavg_accumul_' (cfg.egTrialType)]).([cfg.events{iEv} 'L']),3);
 
  if strcmp(cfg.events{iEv},'preturn')
    taxis  = avgDffSumm.timeAxis_preturn;
    yl     = [-1.5 1];
  elseif strcmp(cfg.events{iEv},'turn')
    taxis  = avgDffSumm.timeAxis;
    yl     = [-.1 4];
  else
    taxis  = avgDffSumm.timeAxis;
    yl     = [-.1 3];
  end

  % single-panel summary, contra and ipsi, towers
  avgAT_contra = (avgAT(:,1:2:nROI,1) + avgAT(:,2:2:end,2))./2;
  avgAT_ipsi   = (avgAT(:,1:2:nROI,2) + avgAT(:,2:2:end,1))./2;
  
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  hold(axs,'on')
  plot([0 0],yl,'--','color',[.7 .7 .7])
  plot([taxis(1) taxis(end)],[0 0],'--','color',[.7 .7 .7])
  for iROI = 1:numel(lbls)
    plot(taxis, avgAT_contra(:,iROI),'-', 'LineWidth', 0.75, 'color', cl(iROI,:));
  end
  xlim([taxis(1) taxis(end)])
  ylim(yl)
  xlabel('Time (s)')
  ylabel('\DeltaF/F (z-score)')
  title(sprintf('Contralateral %s',cfg.events{iEv}));

  
  iPanel = iPanel + 1;
  axs    = fig.panel(iPanel);
  hold(axs,'on')
  plot([0 0],yl,'--','color',[.7 .7 .7])
  plot([taxis(1) taxis(end)],[0 0],'--','color',[.7 .7 .7])
  for iROI = 1:numel(lbls)
    plot(taxis, avgAT_ipsi(:,iROI),'-', 'LineWidth', 0.75, 'color', cl(iROI,:));
  end
  xlim([taxis(1) taxis(end)])
  ylim(yl)
  xlabel('Time (s)')
  ylabel('\DeltaF/F (z-score)')
  title(sprintf('Ipsilateral %s',cfg.events{iEv}));

end

end