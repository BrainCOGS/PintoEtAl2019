function Pinto2019_figS6_WFdynamics_epochComp(analysisFilePath)

% Pinto2019_figS6_WFdynamics_epochComp(analysisFilePath)
% plots avg activity scatter plots comparing different maze positions
% analysisFilePath is path for data analysis files to be loaded

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% refer to Pinto2019_fig3_WFdynamics.m
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.ROIcontraOnly    = true; % to plot just contralateral ROIs 
cfg.refPos           = 25; % x axis for plots
cfg.compPos          = 75:50:275; % x axis for plots
cfg.nBins            = 1; % average over this # bins on either side of positions
cfg.trialType        = 'correctL';
cfg.mazes            = {'visGuide','accumul'};
cfg.mazeLbls         = {'navigate to target','accumulating towers'};

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/avgDffSummary.mat'],'avgDffSumm');

%% Configure figure and panels
layout       = [ 1:5   ...
               ; 6:10  ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% Plot scatters
taxis   = avgDffSumm.spaceBins;
ROIlbl  = avgDffSumm.ROIlbl;

for iMz = 1:numel(cfg.mazes)
  
  if iMz == 1; [cl,ROIlbl] = getDefaultROIcl(ROIlbl); end
  dff         = avgDffSumm.(['ROIavg_' cfg.mazes{iMz} '_' cfg.trialType]); 
  if cfg.ROIcontraOnly 
    side   = cfg.trialType(end);
    idx    = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
    ROIlbl = ROIlbl(idx);
    dff    = dff(:,idx,:);
    cl     = cl(idx,:);
  end

  idx1 = find(taxis == cfg.refPos)-cfg.nBins:find(taxis == cfg.refPos)+cfg.nBins;
  c1   = squeeze(nanmean(dff(idx1,:,:)));
  c1r  = mean(c1,2);
  c1s  = std(c1,0,2)./sqrt(size(c1,2)-1);

  for iPos = 1:numel(cfg.compPos)
    idx2 = find(taxis == cfg.compPos(iPos))-cfg.nBins:find(taxis == cfg.compPos(iPos))+cfg.nBins;
    c2   = squeeze(nanmean(dff(idx2,:,:)));
    c2r  = mean(c2,2);
    c2s  = std(c2,0,2)./sqrt(size(c2,2)-1);
    
    % plot
    iPanel = iPanel + 1;
    axs    = fig.panel(iPanel);
    hold(axs,'on')
    
    for iP = 1:size(c1r,1)
      plot([c1r(iP)-c1s(iP) c1r(iP)+c1s(iP)],[c2r(iP) c2r(iP)], ...
        '-','linewidth',.75,'color',cl(iP,:))
      plot([c1r(iP) c1r(iP)],[c2r(iP)-c2s(iP) c2r(iP)+c2s(iP)], ...
        '-','linewidth',.75,'color',cl(iP,:))
    end
    
    axis tight
    yl     = get(axs,'ylim');
    xl     = get(axs,'xlim');
    xlim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]);
    ylim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]);
    yl     = get(axs,'ylim');
    plot(yl,yl,'--','color',[.7 .7 .7])
    xlabel(['\DeltaF/F (y = ' num2str(cfg.refPos) ')'])
    ylabel(['\DeltaF/F (y = ' num2str(cfg.compPos(iPos)) ')'])

    if iPos == 3; title(cfg.mazeLbls{iMz}); end
  end
end



%% Save figure to disk
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end
