function stats = Pinto2019_figS8_WFdecodingViewAngCtrl(analysisFilePath,summaryFile)

% Pinto2019_figS8_WFdecodingViewAngCtrl(analysisFilePath,summaryFile)
% plots decoding results for widefield data
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To perform decoding for each recording (repo: widefieldImaging) 
%    (note that only 1.2, 1.8 are used in this figure but summary function requires all to run)
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
%    1.5) dffPxlDecoderLeaveOneOut(recpath,spockFlag,decodeWhat,'ridge',0)
%         This will do the ROI-based decoders excluding ROI pairs
%    1.6) dffPxlDecoderLeaveOneOut(recpath,spockFlag,decodeWhat,'ridge',1)
%         This will do the ROI-based decoders excluding ROI pairs, with view angle correction
%    1.7) dffPxlDecoderFromSingleROI(recpath,spockFlag,decodeWhat,'ridge',0)
%         This will do the pxl-based decoders separately for each ROI
%    1.8) dffPxlDecoderFromSingleROI(recpath,spockFlag,decodeWhat,'ridge',1)
%         This will do the pxl-based decoders separately for each ROI with view angle correction
% 2) To summarize experiments: (repo: widefieldImaging)
%    summarizeDecoding_reduced([],[],whichDecoder)
%    whichDecoder is a string that should be set to 'choice', 'prevchoice', and 'evidence'
%    (will save decodingSummary_*.mat, loaded here)
% 3) To compare accuracy across bin sizes: compareDecoderBinning.m (repo: imaging)
%    (will save dffPxlDecoder_binComp.mat, loaded here)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.whichDecoder     = 'allPxls'; 
cfg.decodeAccComp    = {'allPxls_vAngResid','ROIpxls_vAngResid'}; 
cfg.accCompWhat      = {'accuracy_max'}; % ,'accuracy_firsty'
cfg.decWeightGroups  = {'evidence','choice','prevchoice'};%,'viewangle'};
cfg.decWeightLbls    = {'evidence','choice','prev. choice'};%,'view angle'};
cfg.whichStats       = 'allPxls_vAngResid_ROIpxls_vAngResid';

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
load([analysisFilePath '/dffPxlDecoder_binComp.mat'],'decoderComp');
decoderBinComp = decoderComp;

%% cluster
load(summaryFile,'decodeStats');
stats = decodeStats;

%% Configure figure and panels
layout       = [ 1 2 3 4 ...
               ; 5 6 7 8 
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% Evidence decoding accuracy across bins
[iPanel,decoderBinComp]  = pnlBins(fig, iPanel, decoderBinComp, 'evidence');

%% Evidence decoding accuracy across bins
[iPanel,decoderBinComp]  = pnlBins(fig, iPanel, decoderBinComp, 'choice');

%% Evidence decoding accuracy + weights, all ROIs
[iPanel,stats.evidenceDecoder]  = pnlDecodeSummary(fig, iPanel, decodeSumm_evidence, cfg, 'evidence', stats.evidenceDecoder);

%% Choice decoding accuracy + weights, all ROIs
[iPanel,stats.choiceDecoder]  = pnlDecodeSummary(fig, iPanel, decodeSumm_choice, cfg, 'choice', stats.choiceDecoder);

%% Choice decoding accuracy + weights, all ROIs
[iPanel,stats.prevchoiceDecoder]  = pnlDecodeSummary(fig, iPanel, decodeSumm_prevchoice, cfg, 'prevchoice', stats.prevchoiceDecoder);

%% save analysis summary file

if ~isempty(summaryFile)
  decodeStatsViewAngCtrl = stats;
  if isempty(dir(summaryFile))
    save(summaryFile,'decodeStatsViewAngCtrl','decoderBinComp','-v7.3')
  else
    save(summaryFile,'decodeStatsViewAngCtrl','decoderBinComp','-append')
  end
end
%% Save figure to disk
% versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% plot accuracyper bin
function [iPanel,decoderBinComp] = pnlBins(fig, iPanel, decoderBinComp, whichDecoder)

%% plot accuracy vs. space
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

% data
scale   = (1000/widefieldParams.pxlPerMM) * widefieldParams.dsFactor;
xaxis   = decoderBinComp.spatialBinFactor * scale ;
acc     = decoderBinComp.([whichDecoder '_accMax']);
accm    = nanmean(acc);
accs    = nanstd(acc,0,1)./sqrt(size(acc,1)-1);
thiscl  = [0 0 0];

errbar(xaxis,accm,accs,thiscl,.75,0,'-');
plot(xaxis,accm,'.','color',thiscl,'markersize',7);

allacc = acc;

% data roi
xaxis   = mean(decoderBinComp.ROIsize_avg)*scale;
acc     = decoderBinComp.([whichDecoder '_accMax_ROI']);
accm    = nanmean(acc);
accs    = nanstd(acc,0,1)./sqrt(size(acc,1)-1);
thiscl  = [1 0 0];

errbar(xaxis,accm,accs,thiscl,.75,0,'-');
plot(xaxis,accm,'.','color',thiscl,'markersize',7);

% stats
allacc        = [allacc acc];
[nmice,nbins] = size(allacc);
sizevec       = repmat(1:nbins,[nmice 1]);
micevec       = repmat((1:nmice)',[1 nbins]);

[thisp,thistable,thisstats] ...
              = anovan(allacc(:),{sizevec(:),micevec(:)},'varnames',{'pxl size','mice'},'display','off');
thismultcomp  = multcompare(thisstats,'display','off');

decoderBinComp.([whichDecoder '_stats']).ANOVAp        = thisp;
decoderBinComp.([whichDecoder '_stats']).ANOVAtable    = thistable;
decoderBinComp.([whichDecoder '_stats']).ANOVAstats    = thisstats;
decoderBinComp.([whichDecoder '_stats']).ANOVAmultComp = thismultcomp;

% labels
switch whichDecoder 
  case 'evidence'
    titlestr = 'Evidence (Cumul. \Delta)';
    ylbl     = 'Cross-val. acc. (r)';
    yl       = [0 .5];

  case 'choice'
    titlestr = 'Upcoming choice';
    ylbl     = 'Cross-val. acc. (P correct)';
    yl       = [.5 .85];
end

ylabel(ylbl)
xlabel('Pixel size (um)')
title(titlestr)
xlim([250 2500])
ylim(yl)
set(axs,'xscale','log')

text(67,yl(2)*.9,sprintf('p = %1.2g',thisp(1)),'fontsize',9,'color','k')

end

%% ------------------------------------------------------------------------
%% plot accuracy and weights for all-ROI decoder
function [iPanel,stats] = pnlDecodeSummary(fig, iPanel, decodeSumm, cfg, decodeWhat,stats)

%% plot accuracy vs. space
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

% data
xaxis   = decodeSumm.xaxis;
acc     = decodeSumm.(cfg.whichDecoder).accuracy;
accm    = squeeze(nanmean(acc,3));
accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
thiscl  = [0 0 0];

plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)

% data view ang
xaxis   = decodeSumm.xaxis;
acc     = decodeSumm.([cfg.whichDecoder '_vAngResid']).accuracy;
accm    = squeeze(nanmean(acc,3));
accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
thiscl  = [1 0 0];

plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)

% shuffle
acc     = decodeSumm.(cfg.whichDecoder).accuracy_shuffle;
accm    = squeeze(nanmean(acc,3));
accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
thiscl  = [.7 .7 .7];

plot(xaxis,accm,'-','color',thiscl,'linewidth',.75)
plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.25)
plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.25)

% labels
switch decodeWhat
  case 'evidence'
    titlestr = 'Evidence (Cumul. \Delta)';
    ylbl     = 'Cross-val. acc. (r)';
    ylim([-.1 .45])

  case 'choice'
    titlestr = 'Upcoming choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .85])

  case 'prevchoice'
    titlestr = 'Previous choice';
    ylbl     = 'Cross-val. acc. (prop. correct)';
    ylim([.45 .65])
end

ylabel(ylbl)
xlabel('y pos. (cm)')
title(titlestr)
xlim([0 300])

%% Evidence decoding accuracy comp (multi vs single bilateral ROI)
iPanel      = pnlAccComp(fig, iPanel, decodeSumm, cfg, decodeWhat,stats);

end

%% ------------------------------------------------------------------------
%% compare accuracy for different decoders
function iPanel = pnlAccComp(fig, iPanel, decodeSumm, cfg, decodeWhat,stats)

for iComp = 1:numel(cfg.accCompWhat)
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  
  allROIvec    = decodeSumm.(cfg.decodeAccComp{1}).(cfg.accCompWhat{iComp});
  bar(1,nanmean(allROIvec),'facecolor','k','edgecolor','k');
  errorbar(1,nanmean(allROIvec),nanstd(allROIvec)./sqrt(numel(decodeSumm.mice)-1),'k-');
  plotlbls{1}  = 'all ROIs';

  datavec      = decodeSumm.(cfg.decodeAccComp{2}).(cfg.accCompWhat{iComp});
  nROI         = size(datavec,2);
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
  end

  set(axs,'xtick',1:nROI+1,'xticklabel',plotlbls)
  rotateXLabels(axs,90)
  
  % labels
  switch decodeWhat
    case 'evidence'
      titlestr = 'Evidence (Cumul. \Delta)';
      ylbl     = 'cross-val. acc. (r)';
      switch cfg.accCompWhat{iComp}
        case 'accuracy_max'
          ylim([0 .5]);
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
          ylim([.5 .9]);
        case 'accuracy_avg'
          ylim([.5 .55]);
        otherwise
          ylim([0 200]);
      end

    case 'prevchoice'
      titlestr = 'Previous choice';
      ylbl     = 'cross-val. acc. (prop. correct)';
      switch cfg.accCompWhat{iComp}
        case 'accuracy_max'
          ylim([.5 .65]);
        case 'accuracy_avg'
          ylim([.5 .56]);
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
