function stats = Pinto2019_fig6_WFdecoding(analysisFilePath,summaryFile)

% Pinto2019_fig6_WFdecoding(analysisFilePath,summaryFile)
% plots decoding results for widefield data
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To perform decoding for each recording (repo: widefieldImaging) 
%    (note that only 1.1, 1.3 and 1.5 are used in this figure)
%    1.1) dffPxlDecoder(recpath,spockFlag,decodeWhat,'ridge',0,0,[],1)
%         where recpath is full path for recording folder, spockFlag is
%         true if running on spock (ie PNI cluster), decodeWhat is a string
%         that should be set to 'choice', 'prevchoice', and 'evidence'
%         This will do the ROI-based decoders
%    1.2) dffPxlDecoder(recpath,spockFlag,decodeWhat,'ridge',1,0,[],1)
%         This will do the ROI-based decoders with view angle correction
%    1.3) dffPxlDecoder(recpath,spockFlag,decodeWhat,'ridge',0,0,[],0)
%         This will do the pxl-based decoders 
%    1.4) dffPxlDecoder(recpath,spockFlag,decodeWhat,'ridge',1,0,[],0)
%         This will do the pxl-based decoders with view angle correction
%    1.5) dffPxlDecoderFromSingleROI(recpath,spockFlag,decodeWhat,'ridge',0)
%         This will do the pxl-based decoders separately for each ROI
%    1.6) dffPxlDecoderFromSingleROI(recpath,spockFlag,decodeWhat,'ridge',1)
%         This will do the pxl-based decoders separately for each ROI with view angle correction
% 2) To summarize experiments: (repo: widefieldImaging)
%    summarizeDecoding_reduced([],[],whichDecoder)
%    whichDecoder is a string that should be set to 'choice', 'prevchoice', and 'evidence'
%    (will save decodingSummary_*.mat, loaded here)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.whichDecoder     = {'allPxls16x','allPxls'}; 
cfg.decodeAccComp    = {'allPxls16x','ROIpxls'}; 
cfg.accCompWhat      = {'accuracy_max'}; 
cfg.decWeightGroups  = {'evidence','choice','prevchoice'};
cfg.decWeightLbls    = {'evidence','choice','prev. choice'};
cfg.egRec            = 'ai3/20170202/';
cfg.egframes         = 9:4:27;
cfg.whichStats       = 'allPxls_ROIpxls';
cfg.filePath         = analysisFilePath;
cfg.plotSingBarP     = false;
cfg.evidenceCl       = widefieldParams.myblue;
cfg.choiceCl         = widefieldParams.mypurple;

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/decodingSummary_evidence.mat'],'decodeSumm');
decodeSumm_evidence   = decodeSumm;
load([analysisFilePath '/decodingSummary_choice.mat'],'decodeSumm');
decodeSumm_choice     = decodeSumm;
load([analysisFilePath '/decodingSummary_prevchoice.mat'],'decodeSumm');
decodeSumm_prevchoice = decodeSumm;

%% compile stats 
stats.evidenceDecoder   = decodeSumm_evidence.(cfg.decodeAccComp{end}).stats; 
stats.evidenceDecoder   = combinedANOVA(decodeSumm_evidence,stats.evidenceDecoder);
stats.choiceDecoder     = decodeSumm_choice.(cfg.decodeAccComp{end}).stats; 
stats.choiceDecoder     = combinedANOVA(decodeSumm_choice,stats.choiceDecoder);
stats.prevchoiceDecoder = decodeSumm_prevchoice.(cfg.decodeAccComp{end}).stats;
stats.prevchoiceDecoder = combinedANOVA(decodeSumm_prevchoice,stats.prevchoiceDecoder);

%% Configure figure and panels
layout       = [ 1  2  3  4  5  6   ...
               ; 7  7  8  8  9  13  ...
               ; 10 10 11 11 12 13  ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;
fig6stats    = [];

%% Decoding weights example
load(summaryFile,'regMaps');
thism  = mouseAndDateFromFileName(cfg.egRec);
rois   = regMaps.gcamp.outline{strcmpi(regMaps.gcamp.mice,thism)};
[iPanel,fig6stats] = pnlEgWeights(fig, iPanel, decodeSumm_evidence, cfg, rois,fig6stats);

%% weight spatial correlation
[iPanel,fig6stats] = pnlSpatialCorr(fig, iPanel, decodeSumm_evidence, decodeSumm_choice, cfg,fig6stats);

%% Evidence decoding accuracy + weights, all ROIs
[iPanel,stats.evidenceDecoder,fig6stats]  = pnlDecodeSummary(fig, iPanel, decodeSumm_evidence, cfg, 'evidence', stats.evidenceDecoder,fig6stats);

%% Choice decoding accuracy + weights, all ROIs
[iPanel,stats.choiceDecoder,fig6stats]    = pnlDecodeSummary(fig, iPanel, decodeSumm_choice, cfg, 'choice', stats.choiceDecoder,fig6stats);

%% (strip down a bit for saving)
decodeSumm_evidence   = rmMouseFields(decodeSumm_evidence);
decodeSumm_choice     = rmMouseFields(decodeSumm_choice);
decodeSumm_prevchoice = rmMouseFields(decodeSumm_prevchoice);
decodeStats           = stats;

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'decodeStats','decodeSumm_evidence','decodeSumm_choice',...
                     'decodeSumm_prevchoice','fig6stats','-v7.3')
  else
    save(summaryFile,'decodeStats','decodeSumm_evidence','decodeSumm_choice',...
                     'decodeSumm_prevchoice','fig6stats','-append')
  end
end

%% Save figure to disk
% versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% plot example weights for all-pxl decoder
function [iPanel, figstats] = pnlEgWeights(fig, iPanel, decodeSumm, cfg, rois, figstats)

load([cfg.filePath 'brainMask'])
mask = imresize(brainMask,1/16);
load([cfg.filePath cfg.egRec 'reg2ref.mat'],'bregma');
bregma = floor(32- bregma ./ 4);

load([cfg.filePath cfg.egRec 'decoder_evidence_ridge_allTrials_space_binned4x.mat'],'decoder');
egweights = decoder.weights(:,:,cfg.egframes);
cmap      = red2blue;
cmap(1,:) = [0 0 0];
cl        = [-max(abs(egweights(:)))-.0001 max(abs(egweights(:)))+.0001];

for iFrame = 1:numel(cfg.egframes)
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  
  thisf       = egweights(:,:,iFrame);
  thisf(mask) = nan;
  imagesc(flipud(thisf),cl)
  colormap(cmap);
  axis image; axis off; axis xy
  plot(bregma(2)-1,bregma(1)+2,'k+')
  text(1,35,sprintf('y = %d cm',decodeSumm.xaxis(cfg.egframes(iFrame))))
  axis image; axis off; axis xy
  
  for iROI = 1:numel(rois)
    plot(1+((rois{iROI}(:,2))./4),32-((rois{iROI}(:,1))./4),'k-','linewidth',.5);
  end

  if iFrame == numel(cfg.egframes)
    cbar              = smallcolorbar(axs);
    pos               = cbar.Position;
    cbar.Units        = get(axs,'units');
    cbar.Position     = pos;
    cbar.Label.String = 'Dec. weight (a.u.)';
  end
end

figstats.panelA.vasculatureMask      = mask;
figstats.panelA.bregmaCoordinates    = bregma;
figstats.panelA.ROIoutlines          = rois;
figstats.panelA.decodingWeightImages = egweights;
figstats.panelA.decoderYpositions    = decodeSumm.xaxis(cfg.egframes);
figstats.panelA.dataLabel            = 'Evidence decoding weight (a.u.)';

end

%% ------------------------------------------------------------------------
%% plot example weights for all-pxl decoder
function [iPanel, figstats] = pnlSpatialCorr(fig, iPanel, decodeSumm_evidence, decodeSumm_choice, cfg, figstats)

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

nmice = size(decodeSumm_evidence.allPxls.weightSpatialCorr,2);
thisx = decodeSumm_evidence.allPxls.pxlDistAxis;
thism = nanmean(decodeSumm_evidence.allPxls.weightSpatialCorr,2);
thiss = nanstd(decodeSumm_evidence.allPxls.weightSpatialCorr,0,2)./sqrt(nmice-1);
h(1)  = plot(thisx(2:end),thism(2:end),'-','linewidth',1,'color',cfg.evidenceCl); 
plot(thisx(2:end),thism(2:end)-thiss(2:end),'--','linewidth',.5,'color',cfg.evidenceCl); 
plot(thisx(2:end),thism(2:end)+thiss(2:end),'--','linewidth',.5,'color',cfg.evidenceCl); 

figstats.panelB.xaxis              = thisx(2:end);
figstats.panelB.yaxis_evidence_avg = thism(2:end);
figstats.panelB.yaxis_evidence_sem = thiss(2:end);

thisx = decodeSumm_choice.allPxls.pxlDistAxis;
thism = nanmean(decodeSumm_choice.allPxls.weightSpatialCorr,2);
thiss = nanstd(decodeSumm_choice.allPxls.weightSpatialCorr,0,2)./sqrt(nmice-1);
h(2)  = plot(thisx(2:end),thism(2:end),'-','linewidth',1,'color',cfg.choiceCl); 
plot(thisx(2:end),thism(2:end)-thiss(2:end),'--','linewidth',.5,'color',cfg.choiceCl); 
plot(thisx(2:end),thism(2:end)+thiss(2:end),'--','linewidth',.5,'color',cfg.choiceCl); 

figstats.panelB.yaxis_choice_avg = thism(2:end);
figstats.panelB.yaxis_choice_sem = thiss(2:end);

thism = nanmean(decodeSumm_evidence.allPxls.weightSpatialCorr_shuffle,2);
thiss = nanstd(decodeSumm_evidence.allPxls.weightSpatialCorr_shuffle,0,2)./sqrt(nmice-1);
h(3)  = plot(thisx(2:end),thism(2:end),'-','linewidth',1,'color',[.7 .7 .7]); 
plot(thisx(2:end),thism(2:end)-thiss(2:end),'--','linewidth',.5,'color',[.7 .7 .7]); 
plot(thisx(2:end),thism(2:end)+thiss(2:end),'--','linewidth',.5,'color',[.7 .7 .7]); 

figstats.panelB.yaxis_shuffle_avg = thism(2:end);
figstats.panelB.yaxis_shuffle_sem = thiss(2:end);

figstats.panelB.xlabel = 'ROI dist (um)';
figstats.panelB.ylabel = 'Weight correlation (r)';

legend(h,{'evid.','choice','shuff.'},'location','northeast')
legend('boxoff')

xlabel('ROI dist (um)')
ylabel('Weight correlation (r)')

xlim([250 3000])
ylim([-.05 .3])

end

%% ------------------------------------------------------------------------
%% plot accuracy for decoder
function [iPanel,stats, figstats] = pnlDecodeSummary(fig, iPanel, decodeSumm, cfg, decodeWhat,stats, figstats)

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

% data
xaxis   = decodeSumm.xaxis;
acc     = decodeSumm.allROIs.accuracy;
accm    = squeeze(nanmean(acc,3));
accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
thiscl  = [1 0 0];

plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)

switch decodeWhat
  case 'evidence'
    panellbl = 'panelC';
    figstats.(panellbl).ylabel = 'Cross-val. decoding acc. (r, data vs. pred)';
    figstats.(panellbl).title  = 'Evidence (Cumul. \Delta)';
  case 'choice'
    panellbl = 'panelF';
    figstats.(panellbl).ylabel = 'Cross-val. decoding acc. (P correct)';
    figstats.(panellbl).title  = 'Upcoming choice';
end

figstats.(panellbl).xaxis          = xaxis;
figstats.(panellbl).yaxis_ROIs_avg = accm;
figstats.(panellbl).yaxis_ROIs_sem = accs;
    
% data
for iDecoder = 1:numel(cfg.whichDecoder)
  xaxis   = decodeSumm.xaxis;
  acc     = decodeSumm.(cfg.whichDecoder{iDecoder}).accuracy;
  accm    = squeeze(nanmean(acc,3));
  accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
  thiscl  = [0 0 0]+(iDecoder-1)*.5;

  plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
  plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
  plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)
  
  if iDecoder == 1
    figstats.(panellbl).yaxis_pxls1120um_avg = accm;
    figstats.(panellbl).yaxis_pxls1120um_sem = accs;
  else
    figstats.(panellbl).yaxis_pxls280um_avg = accm;
    figstats.(panellbl).yaxis_pxls280um_sem = accs;
  end
end

% shuffle
acc     = decodeSumm.(cfg.whichDecoder{1}).accuracy_shuffle;
accm    = squeeze(nanmean(acc,3));
accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
thiscl  = [.8 .8 .8];

plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)

figstats.(panellbl).yaxis_shuffle_avg = accm;
figstats.(panellbl).yaxis_shuffle_sem = accs;
figstats.(panellbl).xlabel = 'ypos (cm)';

titlestr = 'Evidence (Cumul. \Delta)';
ylbl     = 'Cross-val. decoding acc. (r, data vs. pred)';
ylim([-.1 .5])

ylabel(ylbl)
xlabel('y pos. (cm)')
title(titlestr)
xlim([0 300])

% labels
switch decodeWhat
  case 'evidence'
    titlestr = 'Evidence (Cumul. \Delta)';
    ylbl     = 'Cross-val. acc. (r)';
    ylim([-.1 .5])

  case 'choice'
    titlestr = 'Upcoming choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .8])

  case 'prevchoice'
    titlestr = 'Previous choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .75])
end

ylabel(ylbl)
xlabel('y pos. (cm)')
title(titlestr)
xlim([0 300])

%% Evidence decoding accuracy comp (multi vs single bilateral ROI)
[iPanel,stats, figstats]      = pnlAccComp(fig, iPanel, decodeSumm, cfg, decodeWhat, stats, figstats);

end

%% ------------------------------------------------------------------------
%% compare accuracy for different decoders
function [iPanel,stats, figstats] = pnlAccComp(fig, iPanel, decodeSumm, cfg, decodeWhat, stats, figstats)

%% average time course without error bars
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

xaxis   = decodeSumm.xaxis;
acc  = decodeSumm.(cfg.decodeAccComp{2}).accuracy;
accm = squeeze(nanmean(acc,3));
nROI = size(accm,1);
if nROI == 8; ROIlbl = decodeSumm.ROIlbl_unil;
else; ROIlbl = decodeSumm.ROIlbl; end
for iROI = 1:nROI
  cl = getDefaultROIcl(ROIlbl{iROI});
  plot(xaxis,accm(iROI,:),'-','color',cl,'linewidth',.75)
end
% labels
switch decodeWhat
  case 'evidence'
    titlestr = 'Evidence (Cumul. \Delta)';
    ylbl     = 'Cross-val. acc. (r)';
    ylim([-.1 .5])

  case 'choice'
    titlestr = 'Upcoming choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .8])

  case 'prevchoice'
    titlestr = 'Previous choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .75])
end

ylabel(ylbl)
xlabel('y pos. (cm)')
xlim([0 300])

switch decodeWhat
  case 'evidence'
    panellbl = 'panelD';
    figstats.(panellbl).ylabel = 'Cross-val. decoding acc. (r, data vs. pred)';
    figstats.(panellbl).title  = 'Evidence (Cumul. \Delta)';
  case 'choice'
    panellbl = 'panelG';
    figstats.(panellbl).ylabel = 'Cross-val. decoding acc. (P correct)';
    figstats.(panellbl).title  = 'Upcoming choice';
end

figstats.(panellbl).xaxis           = xaxis;
figstats.(panellbl).yaxis_ROIs      = accm;
figstats.(panellbl).yaxis_ROIlabels = ROIlbl;
figstats.(panellbl).datanote        = 'rows are ROIs columns are y position';


switch decodeWhat
  case 'evidence'
    panellbl = 'panelE';
    figstats.(panellbl).ylabel = 'Max cross-val. decoding acc. (r, data vs. pred)';
    figstats.(panellbl).title  = 'Evidence (Cumul. \Delta)';
  case 'choice'
    panellbl = 'panelH';
    figstats.(panellbl).ylabel = 'Max cross-val. decoding acc. (P correct)';
    figstats.(panellbl).title  = 'Upcoming choice';
end

figstats.(panellbl).xaxis_dataLabels = [{'all areas'} ROIlbl];

for iComp = 1:numel(cfg.accCompWhat)
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  
%   % choice should always be max for clarity of display
%   if strcmpi(decodeWhat,'choice') && ~isempty(strfind(cfg.accCompWhat{iComp},'avg'))
%     cfg.accCompWhat{iComp} = 'accuracy_max';
%   end
  
%   allPxlvec    = decodeSumm.(cfg.decodeAccComp{1}).(cfg.accCompWhat{iComp});
%   bar(1,nanmean(allPxlvec),'facecolor',[.5 .5 .5],'edgecolor',[.5 .5 .5]);
%   errorbar(1,nanmean(allPxlvec),nanstd(allPxlvec)./sqrt(numel(decodeSumm.mice)-1),'-','color',[.5 .5 .5]);
%   plotlbls{1}  = 'all pxls';

  allROIvec    = decodeSumm.(cfg.decodeAccComp{1}).(cfg.accCompWhat{iComp});
  bar(1,nanmean(allROIvec),'facecolor','k','edgecolor','k');
  errorbar(1,nanmean(allROIvec),nanstd(allROIvec)./sqrt(numel(decodeSumm.mice)-1),'k-');
  plotlbls{1}  = 'all areas';
  
  figstats.(panellbl).yaxis_avg(1) = nanmean(allROIvec);
  figstats.(panellbl).yaxis_sem(1) = nanstd(allROIvec)./sqrt(numel(decodeSumm.mice)-1);
  
  %% stats
  if strcmpi(decodeWhat,'evidence')
    th       = 0;
  else
    th       = .5;
  end
  [~,thisp]    = ttest(allROIvec,th);
  datavec      = decodeSumm.(cfg.decodeAccComp{2}).(cfg.accCompWhat{iComp});
  nROI         = size(datavec,2);
  
  ROIps      = zeros(1,nROI);
  for iROI = 1:nROI
    [~,ROIps(iROI)] = ttest(datavec(:,iROI),th);
  end
  allps      = [thisp ROIps];    
  isSig      = FDR(allps');
  stats.PvalVSchance       = allps;
  stats.PvalVSchance_isSig = isSig;
  
  %%
  if cfg.plotSingBarP 
    if thisp >= .05 || ~isSig(1); txtstr = 'n.s.'; end
    if thisp < .05 && isSig(1);   txtstr = '*';    end
    if thisp < .01 && isSig(1);   txtstr = '**';   end
    if thisp < .001 && isSig(1);  txtstr = '***';  end
    thisy = nanmean(allROIvec)+nanstd(allROIvec)./sqrt(numel(decodeSumm.mice)-1)+.05;
    text(1,thisy,txtstr,'color','k','horizontalAlignment','center')
  end
  
  if nROI == 8
    roilbls = decodeSumm.ROIlbl_unil;
  else
    roilbls = decodeSumm.ROIlbl;
  end
  for iROI = 1:nROI
    thismean = mean(datavec(:,iROI));
    thissem  = std(datavec(:,iROI))./sqrt(numel(decodeSumm.mice)-1);
    [cl,plotlbls{end+1}] = getDefaultROIcl(roilbls{iROI});
    bar(iROI+1,thismean,'edgecolor',cl,'facecolor',cl);
    errorbar(iROI+1,thismean,thissem,'-','color',cl);
    
    if cfg.plotSingBarP
      if allps(iROI+1) >= .05 || ~isSig(iROI+1); txtstr = 'n.s.'; end
      if allps(iROI+1) < .05 && isSig(iROI+1);   txtstr = '*';    end
      if allps(iROI+1) < .01 && isSig(iROI+1);   txtstr = '**';   end
      if allps(iROI+1) < .001 && isSig(iROI+1);  txtstr = '***';  end
      thisy = thismean+thissem+.05;
      text(iROI+1,thisy,txtstr,'color','k','horizontalAlignment','center')
    end
    
    figstats.(panellbl).yaxis_avg(iROI+1) = thismean;
    figstats.(panellbl).yaxis_sem(iROI+1) = thissem;
  end
  stats.dataLbls = plotlbls;
  set(axs,'xtick',1:nROI+1,'xticklabel',plotlbls)
  rotateXLabels(axs,90)
  
  % labels
  switch decodeWhat
    case 'evidence'
      titlestr = 'Evidence (Cumul. \Delta)';
      ylbl     = 'cross-val. acc. (r)';
      switch cfg.accCompWhat{iComp}
        case 'accuracy_max'
          ylim([0 .6]);
        case 'accuracy_avg'
          ylim([0 .2]);
        otherwise
          ylim([0 200]);
      end
    case 'choice'
      titlestr = 'Upcoming choice';
      ylbl     = 'cross-val. acc. (prop. correct)';
      switch cfg.accCompWhat{iComp}
        case 'accuracy_max'
          ylim([.5 .8]);
        case 'accuracy_avg'
          ylim([.5 .56]);
        otherwise
          ylim([0 200]);
      end

    case 'prevchoice'
      titlestr = 'Previous choice';
      ylbl     = 'cross-val. acc. (prop. correct)';
      switch cfg.accCompWhat{iComp}
        case 'accuracy_max'
          ylim([.5 .7]);
        case 'accuracy_avg'
          ylim([.5 .6]);
        otherwise
          ylim([0 200]);
      end
  end

  switch cfg.accCompWhat{iComp}
    case 'accuracy_max'
      ylbl = sprintf('%s\nMax %s',titlestr,ylbl);
    case 'accuracy_avg'
      ylbl = sprintf('%s\nAvg %s',titlestr,ylbl);
    case 'accuracy_firsty'
      ylbl = sprintf('%s\nDecode latency (cm)',titlestr);
    case 'accuracy_maxy'
      ylbl = sprintf('%s\nPos. max acc. (cm)',titlestr);
  end
  ylabel(ylbl);
  
  yl   = get(axs,'ylim');
  pval = stats.(cfg.whichStats).(cfg.accCompWhat{iComp}).anovaP(1);
  plot([1 nROI],[yl(2) yl(2)],'k-')
  text(nROI/2+1,yl(2)*1.02,sprintf('P = %1.2g',pval),'color','k','fontsize',8,'horizontalAlignment','center')
end

end

%% strip down data structures for saving
function decodeSumm = rmMouseFields(decodeSumm)

fls = fieldnames(decodeSumm);
for iF = 1:numel(fls)
  if isfield(decodeSumm.(fls{iF}),'mouse')
    decodeSumm.(fls{iF}) = rmfield(decodeSumm.(fls{iF}),'mouse');
  end
end
end

%% ANOVA combining decoders
function stats = combinedANOVA(decodeSumm,stats)
anova_types    = {'accuracy_avg','accuracy_max','accuracy_maxy','accuracy_firsty'};
decoderComb{1} = {'allPxls','ROIpxls'};
decoderComb{2} = {'allPxls_vAngResid','ROIpxls_vAngResid'};

for iComb = 1:numel(decoderComb)
  for iType = 1:numel(anova_types)
    roivec       = [];
    counter      = 1;
    name         = [];
    for iDec = 1:numel(decoderComb{iComb})
      name         = [name '_' decoderComb{iComb}{iDec}];
      datavec      = decodeSumm.(decoderComb{iComb}{iDec}).(anova_types{iType});
      [nmice,nROI] = size(datavec);
      datavec      = datavec(:);
      micevec      = repmat((1:nmice)',[nROI 1]);

      for iROI = 1:nROI
        roivec     = [roivec; ones(nmice,1)*counter];
        counter    = counter+1;
      end
    end
    name = name(2:end);
    [stats.(name).(anova_types{iType}).anovaP,     ...
     stats.(name).(anova_types{iType}).anovaTable, ...
     stats.(name).(anova_types{iType}).anovaStats] ...
                = anovan(datavec(~isnan(datavec)),{roivec(~isnan(datavec)),micevec(~isnan(datavec))},'display','off');
    stats.(name).(anova_types{iType}).multComp     ...
                 = multcompare(stats.(name).(anova_types{iType}).anovaStats,'display','off');
  end
end 

end