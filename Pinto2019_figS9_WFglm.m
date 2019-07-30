function Pinto2019_figS9_WFglm(analysisFilePath,summaryFile)

% Pinto2019_figS9_WFglm(analysisFilePath,summaryFile)
% plots GLM results for widefield data
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To run GLM for each recording session: (repo: widefieldImaging) 
%    dffGLM(recpath), where recpath is full path for recording folder
% 2) To summarize experiments:
%    summarizeGLM.m (will save glmSummary_ROI.mat, loaded here)(repo: widefieldImaging)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.GLMcontraOnly    = true;
cfg.clustMaxPC       = 3; % number of PCs to use for clustering
cfg.clustMaxK        = 5; % maximal number of clusters to test
cfg.clustCl          = [60 179 113; 150 110 70; 10 35 140; 230 186 0]./255;
cfg.clustOrder       = [2 1 4 3];
cfg.useCorrClust     = true;
cfg.plotGLMWeights   = {'tow_R','ch','y'};
cfg.plotWeightLbls   = {'R tower','Choice','y pos'};
cfg.GLMweightGroups  = {'towers','\Delta towers','y pos','kinematic','choice','prev. choice','prev. rw'};
cfg.GLMpredEg        = 'ai3/20170202';
cfg.GLMpredEgNTrials = 8;
cfg.GLMpredEgROI     = {'VISp-L','VISa-L','RSP-L','mMOs-L','MOp-L'};
cfg.GLMfn            = 'dffGLM_space_ridge_ROI.mat';
cfg.egTrialType      = 'R';
cfg.corrCMap         = gray;

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/glmSummary_ROI.mat'],'glmSumm');

%% compile stats 
stats.GLM     = glmSumm.stats;  

%% use corr cluster?
if cfg.useCorrClust
  load(summaryFile,'WFcorrStats');
  stats.clusterGLM = WFcorrStats.clusterCorr;
  % reorder cluster numbering for plotting
  if ~isempty(cfg.clustOrder)
    kIDold = stats.clusterGLM.clustID;
    for iK = 1:stats.clusterGLM.nK
      stats.clusterGLM.clustID(kIDold == cfg.clustOrder(iK)) = iK;
    end
  end
else
  stats.clusterGLM = [];
end

%% retrieve egs
egGLM         = load([analysisFilePath cfg.GLMpredEg '/' cfg.GLMfn]);

%% Configure figure and panels
layout       = [ 1  1  2  ...
               ; 3  4  5  ...
               ; 6  6  6  ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;


%% GLM: prediction eg
iPanel          = pnlGLMPredEg(fig, iPanel, egGLM, cfg);

%% GLM: prediction quantification
[iPanel,stats]  = pnlGLMPredAvg(fig, iPanel, glmSumm, cfg, stats);

%% GLM: weights
[iPanel,stats]  = pnlGLMWeights(fig, iPanel, glmSumm, cfg, stats);

%% GLM weights by cluster
[iPanel,stats]  = pnlGLMWeightByClust(fig, iPanel, glmSumm, cfg, stats);

%% (strip down a bit for saving)
glmSumm       = rmfield(glmSumm,'mouse');
glmStats      = stats;

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'glmSumm','glmStats','-v7.3')
  else
    save(summaryFile,'glmSumm','glmStats','-append')
  end
end
%% Save figure to disk
versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end


%% plot GLM g.o.f.
function [iPanel,stats] = pnlGLMPredAvg(fig, iPanel, glmSumm, cfg, stats)

acc_avg     = glmSumm.stats.accuracy_avg;
acc_sem     = glmSumm.stats.accuracy_sem;
shuffleMean = nanmean(glmSumm.stats.accuracy_shuffle_avg);
ROIlbl      = glmSumm.ROIlbl;
[cl,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.GLMcontraOnly 
  side    = cfg.egTrialType(end);
  idx     = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl  = ROIlbl(idx);
  cl      = cl(idx,:);
  acc_avg = acc_avg(idx);
  acc_sem = acc_sem(idx);
end


iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

for iROI = 1:numel(ROIlbl)
  bar(iROI,acc_avg(iROI),'edgecolor',cl(iROI,:),'facecolor',cl(iROI,:))
  errorbar(iROI,acc_avg(iROI),acc_sem(iROI),'-','color',cl(iROI,:))
end

plot([0 numel(ROIlbl)+1],[shuffleMean shuffleMean])

set(axs,'xtick',1:numel(ROIlbl),'xticklabel',ROIlbl)
rotateXLabels(axs,90)
ylabel('Cross-val. prediction r')

ylim([0 .4])

%% ANOVA
acc          = glmSumm.accuracy';
[nROI,nmice] = size(acc);
acc          = acc(:);
mid          = ones(nROI,nmice)*diag(1:nmice);
mid          = mid(:);
roiid        = repmat((1:nROI)',[nmice 1]);

[stats.GLMaccByROI.anovaP,stats.GLMaccByROI.anovaTable,stats.GLMaccByROI.anovaStats] ...
                           = anovan(acc,{roiid,mid},'display','off');
stats.GLMaccByROI.multComp = multcompare(stats.GLMaccByROI.anovaStats,'display','off');

yl = get(axs,'ylim');
xl = get(axs,'xlim');
text(diff(xl)/2,yl(2)*.98,sprintf('P = %1.2g',stats.GLMaccByROI.anovaP(1)),'horizontalAlignment','center')

    
end

%% plot GLM prediction vs data example
function iPanel = pnlGLMPredEg(fig, iPanel, glmEg, cfg)

% pick 5 trials with best avg prediction for desired ROIs
ROIlbl  = glmEg.dffFit.ROIlbl;
ROIidx  = cellfun(@(x)(find(strcmpi(ROIlbl,x))),cfg.GLMpredEgROI);
[~,lbl] = getDefaultROIcl(cfg.GLMpredEgROI);
bp_y    = cell2mat(glmEg.dffFit.bestpred_y);
bp_yh   = cell2mat(glmEg.dffFit.bestpred_yhat);
data    = bp_y(:,ROIidx);
pred    = bp_yh(:,ROIidx);
ypos    = glmEg.dffFit.cfg.posBins;

nt      = numel(ypos)-1;
nt_eg   = cfg.GLMpredEgNTrials*nt;
tcorr   = zeros(cfg.GLMpredEgNTrials,1);
for iW = 1:nt_eg:size(data,1)
  thisy       = data(iW:min([iW+nt_eg-1 size(data,1)]),:);
  thisyhat    = pred(iW:min([iW+nt_eg-1 size(data,1)]),:);
  tcorr(floor(iW/nt_eg)+1) = corr(thisy(:),thisyhat(:));
end

tidx    = (find(tcorr==max(tcorr),1,'first')-1)*nt:find(tcorr==max(tcorr),1,'first')*nt+nt_eg;
data    = data(tidx,:);
pred    = pred(tidx,:);

% insert nan's between trials for plotting
numNan  = 20;
newData = nan(size(data,1)+numNan*cfg.GLMpredEgNTrials-1,size(data,2));
newPred = newData;
newData(1:nt,:) = data(1:nt,:);
newPred(1:nt,:) = pred(1:nt,:);
for iTrial = 2:cfg.GLMpredEgNTrials
  newData(numNan*(iTrial-1)+nt*(iTrial-1)+1:numNan*(iTrial-1)+nt*iTrial,:) = data(nt*(iTrial-1)+1:nt*iTrial,:);
  newPred(numNan*(iTrial-1)+nt*(iTrial-1)+1:numNan*(iTrial-1)+nt*iTrial,:) = pred(nt*(iTrial-1)+1:nt*iTrial,:);
end
newData(end-numNan+1:end,:) = [];
newPred(end-numNan+1:end,:) = [];

%% plot
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

offset = 3;
for iROI = 1:numel(lbl)
  plot(newData(:,iROI)+(iROI-1)*offset,'k-','linewidth',.5)
  plot(newPred(:,iROI)+(iROI-1)*offset,'-','color',[.5 .5 .5],'linewidth',1)
  text(50,offset*(iROI-1)+1,lbl{iROI})
end

% time scale, dff scale, behavioral variable scales
axis tight; 
axis off
yl = get(axs,'ylim');
plot([-15 100-15],[yl(1)-.01 yl(1)-.01],'k-')
text(50-15, yl(1)-.3, '1 m', 'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8)
plot([-15 -15],[yl(1) yl(1)+offset/3],'k-')
th = text(-45, yl(1)+offset/3/2, '\DeltaF/F (Z-score = 1)', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8);
set(th,'rotation',90)

text(size(newData,1)-200,yl(2)*.98,'Data','fontsize',9,'color','k')
text(size(newData,1)-200,yl(2)*.93,'GLM prediction','fontsize',9,'color',[.5 .5 .5])

end

%% plot GLM weights
function [iPanel, stats] = pnlGLMWeights(fig, iPanel, glmSumm, cfg, stats)

weights     = glmSumm.weights;
nmice       = size(weights,3);
weights_avg = glmSumm.stats.weight_avg;
weights_sem = glmSumm.stats.weight_sem;
predLbls    = glmSumm.predLbls;
ROIlbl      = glmSumm.ROIlbl;
[cl,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.GLMcontraOnly 
  side        = cfg.egTrialType(end);
  idx         = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl      = ROIlbl(idx);
  cl          = cl(idx,:);
  weights_avg = weights_avg(idx,:);
  weights_sem = weights_sem(idx,:);
  weights     = weights(idx,:,:);
end

for iPnl = 1:numel(cfg.plotGLMWeights)
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  
  isPred = arrayfun(@(x)(~isempty(strmatch(cfg.plotGLMWeights{iPnl},x))),predLbls);
  nlags  = sum(isPred);
  iPred  = strcmpi(glmSumm.glmCfg.predList,cfg.plotGLMWeights{iPnl});

  if nlags > 1
    switch glmSumm.cfg.timeOrSpace
      case 'time'
        lags = linspace(-glmSumm.glmCfg.predLagSec{iPred}(1),glmSumm.glmCfg.predLagSec{iPred}(2),nlags);
        xlbl = 'Lag (s)';
      case 'space'
        lags = linspace(-glmSumm.glmCfg.predLagCm{iPred}(1),glmSumm.glmCfg.predLagCm{iPred}(2),nlags);
        xlbl = 'Lag (cm)';
    end
    thisw = [];
    mid   = [];
    roiid = [];
    for iROI = numel(ROIlbl):-1:1
      
      plot(lags,weights_avg(iROI,isPred),'-','linewidth',.75,'color',cl(iROI,:));
      thisw = [thisw; squeeze(sum(weights(iROI,isPred,:),2))];
      mid   = [mid; (1:nmice)'];
      roiid = [roiid; ones(nmice,1)*iROI];
    end
    xlabel(xlbl)
    
    [stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaP,stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaTable,stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaStats] ...
                          = anovan(thisw,{roiid,mid},'display','off');
    stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).multComp = multcompare(stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaStats,'display','off');
    
    yl = get(axs,'ylim');
    xl = get(axs,'xlim');
    text(diff(xl)/2,yl(2)*.98,sprintf('P = %1.2g',stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaP(1)),'horizontalAlignment','center')
  else
    thisw = [];
    mid   = [];
    roiid = [];
    for iROI = 1:numel(ROIlbl)
      bar(iROI,weights_avg(iROI,isPred),'edgecolor',cl(iROI,:),'facecolor',cl(iROI,:))
      errorbar(iROI,weights_avg(iROI,isPred),weights_sem(iROI,isPred),'-','color',cl(iROI,:))
      thisw = [thisw; squeeze(weights(iROI,isPred,:))];
      mid   = [mid; (1:nmice)'];
      roiid = [roiid; ones(nmice,1)*iROI];
    end

    set(axs,'xtick',1:numel(ROIlbl),'xticklabel',ROIlbl)
    rotateXLabels(axs,90) 
    
    [stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaP,stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaTable,stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaStats] ...
                          = anovan(thisw,{roiid,mid},'display','off');
    stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).multComp = multcompare(stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaStats,'display','off');
    
    yl = get(axs,'ylim');
    xl = get(axs,'xlim');
    text(diff(xl)/2,yl(2)*.98,sprintf('P = %1.2g',stats.GLMweightsByROI.(cfg.plotGLMWeights{iPnl}).anovaP(1)),'horizontalAlignment','center')
  end
  
  ylabel(sprintf('%s weight',cfg.plotWeightLbls{iPnl}))
%   axis tight
end

end

%% ------------------------------------------------------------------------
%% average glm weights by cluster (with stats across mice)
function [iPanel,stats] = pnlGLMWeightByClust(fig, iPanel, glmSumm, cfg, stats)

% decide which clustering to use
kID = stats.clusterGLM.clustID;
stats.nK = numel(unique(kID));

% reorder cluster numbering for plotting
kIDold = kID;
for iK = 1:stats.nK
  kID(kIDold == cfg.clustOrder(iK)) = iK;
end

% group by weight types (ie sum(towers), sum(delta), y, avg(kinematic),
% choice, prev choice, etc) 
weights     = glmSumm.weights;
nmice       = size(weights,3);
predLbls    = glmSumm.predLbls;
clustW_all  = cell(numel(cfg.GLMweightGroups),numel(unique(kID)));
clustW_mean = zeros(numel(cfg.GLMweightGroups),numel(unique(kID)));
clustW_sem  = zeros(numel(cfg.GLMweightGroups),numel(unique(kID)));
    
for iPt = 1:numel(cfg.GLMweightGroups)
  switch cfg.GLMweightGroups{iPt}
    case 'towers'
      preds = {'tow_R','tow_L'};
    case {'evidence','\Delta','\Delta towers'}
      preds = {'\Delta'};
    case 'y pos'
      preds = {'y'};
    case 'kinematic'
      preds = {'\theta','d\theta/dt','speed'};
    case 'choice'
      preds = {'ch'};
    case 'prev. choice'
      preds = {'prevch'};
    case 'prev. rw'
      preds = {'prevrw'};
  end
  
  for iK = 1:numel(unique(kID))
    isROI = kID == iK;
    thisw = zeros(sum(isROI),numel(preds),nmice);
    for iPred = 1:numel(preds)
      isPred = arrayfun(@(x)(~isempty(strmatch(preds{iPred},x))),predLbls);
      thisw(:,iPred,:) = sum(weights(isROI,isPred,:),2);
    end
    thisw               = mean(abs(thisw),2);
    clustW_all{iPt,iK}  = squeeze(thisw);
    clustW_mean(iPt,iK) = mean(squeeze(mean(thisw,3)));
    clustW_sem(iPt,iK)  = std(mean(thisw,1),0,3)./sqrt(nmice-1);
  end
end

%% do stats to compare significance of weight differences between clusters
% one-way ANOVA separately for each weight group
stats.GLMweightsByClust.labels = cfg.GLMweightGroups;
for iPred = 1:numel(cfg.GLMweightGroups)
  datavec = [];
  groupid = [];
  mouseid = [];
  for iK = 1:numel(unique(kID))
    datavec = [datavec; mean(clustW_all{iPred,iK})';];
    groupid = [groupid; ones(nmice,1).*iK];
    mouseid = [mouseid; (1:nmice)'];
  end
  [stats.GLMweightsByClust(iPred).anovaP,stats.GLMweightsByClust(iPred).anovaTable,stats.GLMweightsByClust(iPred).anovaStats] ...
                          = anovan(datavec,{groupid,mouseid},'display','off');
  stats.GLMweightsByClust(iPred).multComp = multcompare(stats.GLMweightsByClust(iPred).anovaStats,'display','off');
end

%% plot
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

xt = zeros(1,numel(cfg.GLMweightGroups));
for iPt = 1:numel(cfg.GLMweightGroups)
  xt(iPt) = (iPt-1)*stats.nK + ceil(size(clustW_mean,2)/2) + iPt;
  for iK = 1:numel(unique(kID))
    x = (iPt-1)*stats.nK + iK + (iPt-1);
    if iPt == 1
      lh(iK)     = bar(x,clustW_mean(iPt,iK),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:));
      lh_lbl{iK} = ['clust ' num2str(iK)];
    else
      bar(x,clustW_mean(iPt,iK),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:))
    end
    errorbar(x,clustW_mean(iPt,iK),clustW_sem(iPt,iK),'-','color',cfg.clustCl(iK,:))
  end
  yl = get(axs,'ylim');
  text(xt(iPt),yl(2)*.98,sprintf('P = %1.3g',stats.GLMweightsByClust(iPt).anovaP(1)),'horizontalAlignment','center') % print p val
end

legend(lh,lh_lbl,'location','westoutside'); legend('boxoff')
set(axs, 'xtick', xt, 'xticklabel', cfg.GLMweightGroups)
xlim([0 xt(end)+2])
ylim([0 .5])
ylabel(axs, 'Avg abs. weight (a.u.)')
rotateXLabels(axs,30)
end

