function Pinto2019_fig5_evidenceTuning(analysisFilePath,summaryFile,loadData)

% Pinto2019_fig5_evidenceTuning(analysisFilePath,summaryFile)
% plots evidence breakdown of activity
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) see owf_WFdynamics.m for preprocessing
% 2) To generate per-evidence activity averages: (repo: widefieldImaging)
%    avgDffROIvsEvidence.m (use default inputs)
% 3) To summarize per-evidence activity averages: (repo: widefieldImaging)
%    summarizeAvgDffvsEvidence.m (will save avgDffSummary.mat, loaded here)
%% ------------------------------------------------------------------------

if nargin < 3; loadData = 0; end

%% Analysis configuration
cfg.egAreaRamps1   = 'V1';
cfg.egAreaRamps2   = 'mM2';
cfg.doNormRamps    = false;
cfg.nBootIter      = 50;
cfg.nSDsig         = 2;
cfg.posAxisLim     = [0 200];
cfg.mazeSections   = {[0 200]};%,[200 280],[280 300]};
cfg.fitTasks       = {'accumul','visGuide'};
cfg.fitTrialTypes  = {'correct','error'};
cfg.egAreaFits     = {'V1','V1','mM2'};
cfg.egFitTrialType = {'correct','error','correct'};
cfg.egFitTask      = {'accumul','accumul','accumul'};
cfg.egFitMzSect    = {[0 200],[0 200],[0 200]};
cfg.contraCl       = widefieldParams.myblue;
cfg.contraSh       = widefieldParams.myblueSh;
cfg.ipsiCl         = widefieldParams.myred;
cfg.ipsiSh         = widefieldParams.myredSh;

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});

%% load data
load([analysisFilePath 'avgDffSummaryVsEvidence.mat'],'avgDffSumm')
avgDffvsEvidence = avgDffSumm; clear avgDffSumm;
avgDffvsEvidence = rmfield(avgDffvsEvidence,'mouse');

%% linear fits to evidence curves
if loadData
  load(summaryFile,'avgDffvsEvidence');
else
  avgDffvsEvidence = linEvidFits(avgDffvsEvidence,cfg); 
end

% load evidTunSummaryGP_ROI.mat evidTunSumm
% evidTunSummGP.towers = evidTunSumm; clear evidTunSumm
% evidTunSummGP.towers = rmfield(evidTunSummGP.towers,'mouse');
% 
% load evidTunSummaryGP_ROI_visGuide.mat evidTunSumm
% evidTunSummGP.visGuide = evidTunSumm; clear evidTunSumm
% evidTunSummGP.visGuide = rmfield(evidTunSummGP.visGuide,'mouse');

%% set up figure
layout      = [1 1 4 7 8 9 ...
             ; 2 2 5 7 8 9 ...
             ; 3 3 6 7 8 9 ...
             ];
fig         = PaneledFigure(layout, 'tiny');
iPanel      = 0;

%% example "ramps"
fig5stats=[];
[iPanel,fig5stats]      = pnlROIeg(fig, iPanel, avgDffvsEvidence, cfg.egAreaRamps1, 'accumul', 'correct', false, cfg,fig5stats);
[iPanel,fig5stats]      = pnlROIeg(fig, iPanel, avgDffvsEvidence, cfg.egAreaRamps1, 'accumul', 'error', false, cfg,fig5stats);
[iPanel,fig5stats]      = pnlROIeg(fig, iPanel, avgDffvsEvidence, cfg.egAreaRamps2, 'accumul', 'correct', false, cfg,fig5stats);

%% example lin fits
for iSect = 1:numel(cfg.egFitMzSect)
  [iPanel,fig5stats]    = pnlFits(fig, iPanel, avgDffvsEvidence, cfg, iSect,fig5stats);
end

%% slopes summary
for iSect = 1:numel(cfg.mazeSections)
  [iPanel,fig5stats]    = pnlFitSlopeSummary(fig, iPanel, avgDffvsEvidence, cfg, 'accumul', 'correct', cfg.egFitMzSect{iSect},iSect==1,fig5stats);
end
[iPanel,fig5stats]      = pnlFitSlopeSummary(fig, iPanel, avgDffvsEvidence, cfg, 'accumul', 'error', cfg.egFitMzSect{1},false,fig5stats);
[iPanel,fig5stats]      = pnlFitSlopeSummary(fig, iPanel, avgDffvsEvidence, cfg, 'visGuide', 'correct', cfg.egFitMzSect{1},false,fig5stats);

% %% tuning for either task (old version, gaussian-procees fitting)
% [iPanel,evidTunSummGP.towers]   = pnlGPtuning(fig, iPanel, evidTunSummGP.towers, cfg, 'Accum. towers');
% [iPanel,evidTunSummGP.visGuide] = pnlGPtuning(fig, iPanel, evidTunSummGP.visGuide, cfg, 'Visually guided');
% 
% %% task comparison
% [iPanel,evidTunSummGP] = pnlTaskComp(fig, iPanel, evidTunSummGP);

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'avgDffvsEvidence','fig5stats','-v7.3')
  else
    save(summaryFile,'avgDffvsEvidence','fig5stats','-append')
  end
end

%% Save figure to disk
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ----------------------------------------------
%% plot avg Dff per evidence
function [iPanel,figstats] = pnlROIeg(fig, iPanel, avgDffvsEvidence, egAreaRamps, task, trialType, plotLegend, cfg,figstats)

% get matrices
[~,lbls] = getDefaultROIcl(avgDffvsEvidence.ROIlbl);
isEgROI  = contains(lbls,egAreaRamps);
egEvid   = avgDffvsEvidence.(sprintf('ROIavg_%s_%s',task,trialType))(:,isEgROI,:,:);

% make ipsi negative and contra positive
Ridx               = contains(lbls(isEgROI),'-R');
egEvid(:,Ridx,:,:) = flipud(egEvid(:,Ridx,end:-1:1,:));

% norm between 0 and 1?
if cfg.doNormRamps
  nmice = size(egEvid,4);
  nROI  = size(egEvid,2);
  for iMouse = 1:nmice
    for iROI = 1:nROI
      thismat  = egEvid(:,iROI,:,iMouse);
      thismat  = thismat - min(thismat(:));
      thismat  = thismat ./ max(thismat(:));
      egEvid(:,iROI,:,iMouse) = thismat;
    end
  end
end

% hemisphere, then mouse avg
egEvid  = squeeze(nanmean(nanmean(egEvid,2),4));

% plot
nEv     = size(egEvid,2);
cl      = red2blue(nEv);
iPanel  = iPanel + 1;
axs     = fig.panel(iPanel);
hold(axs,'on')

if strcmp(task,'accumul')
  yl = [-.15 .25];
else
  yl = [-.25 .5];
end

for iSection = 1:numel(cfg.mazeSections)
  plot([cfg.mazeSections{iSection}(1) cfg.mazeSections{iSection}(1)],yl,'--','color',[.7 .7 .7],'linewidth',.25)
  plot([cfg.mazeSections{iSection}(2) cfg.mazeSections{iSection}(2)],yl,'--','color',[.7 .7 .7],'linewidth',.25)
end
for iEv = 1:nEv
  plot(avgDffvsEvidence.spaceBins,egEvid(:,iEv),'-','color',cl(iEv,:),'linewidth',.75)
end

xlim(cfg.posAxisLim);
ylim(yl)

xlabel('y pos (cm)')
ylabel('\DeltaF/F (z-score)')
title(sprintf('%s, %s trials',egAreaRamps,trialType))

if plotLegend
  legend(cellfun(@num2str,num2cell(avgDffvsEvidence.RminusLvals),'UniformOutput',false),'location','best')
end

idxt                                    = find(avgDffvsEvidence.spaceBins==cfg.posAxisLim(1)):find(avgDffvsEvidence.spaceBins==cfg.posAxisLim(end));
figstats.panelA_left.evidenceBinCenters = avgDffvsEvidence.RminusLvals;
figstats.panelA_left.xaxis              = avgDffvsEvidence.spaceBins(idxt);
thislbl                                 = ['data_' egAreaRamps '_' trialType 'Trials_' task];
figstats.panelA_left.(thislbl)          = egEvid(idxt,:);
figstats.panelA_left.xlabel             = 'y pos (cm)';
figstats.panelA_left.ylabel             = '\DeltaF/F (z-score)';
figstats.panelA_left.datanote1          = 'rows are time points and cols are evidence values in data matrix';
figstats.panelA_left.datanote2          = 'negative evidence values are ipsilateral, positive are contralateral';

end

%% ----------------------------------------------
%% fit lines to contra and ipsi evidence, bootstrap for error & significance
function avgDffvsEvidence = linEvidFits(avgDffvsEvidence,cfg)

%% initialize
lbls      = avgDffvsEvidence.ROIlbl;
lblsUnil  = unique(cellfun(@(x)(x(1:end-2)),lbls,'uniformoutput',false),'stable');
nROI      = numel(lblsUnil);
xaxis     = avgDffvsEvidence.RminusLvals';
contraidx = numel(xaxis)/2+1:numel(xaxis);
ipsiidx   = 1:numel(xaxis)/2;

avgDffvsEvidence.linFits.mazeSections = cfg.mazeSections; 
[~,avgDffvsEvidence.linFits.ROIlbl]   = getDefaultROIcl(lblsUnil);
  
%% loop over tasks, trial types, ROIs, maze sections
for iTask = 1:numel(cfg.fitTasks)
  task = cfg.fitTasks{iTask};
  for iTrialType = 1:numel(cfg.fitTrialTypes)
    trialType = cfg.fitTrialTypes{iTrialType};
    if strcmpi(task,'visGuide') && strcmpi(trialType,'error'); continue; end
    for iSection = 1:numel(cfg.mazeSections)
      pos       = cfg.mazeSections{iSection};
      posidx    = avgDffvsEvidence.spaceBins > pos(1) & avgDffvsEvidence.spaceBins <= pos(2);
      for iROI = 1:nROI
        isEgROI  = contains(lbls,lblsUnil{iROI});
        egEvid   = avgDffvsEvidence.(sprintf('ROIavg_%s_%s',task,trialType))(posidx,isEgROI,:,:);
        
        % make ipsi negative and contra positive
        Ridx               = contains(lbls(isEgROI),'-R');
        egEvid(:,Ridx,:,:) = egEvid(:,Ridx,end:-1:1,:);
        nmice              = size(egEvid,4);
        
        % boostrap mice
        contrap1 = zeros(cfg.nBootIter,1);
        ipsip1   = zeros(cfg.nBootIter,1);
        contrap2 = zeros(cfg.nBootIter,1);
        ipsip2   = zeros(cfg.nBootIter,1);
        
        for iBoot = 1:cfg.nBootIter
          midx          = randsample(nmice,nmice,true);
          egEvid_contra = squeeze(nanmean(nanmean(egEvid(:,:,contraidx,midx),2),1));
          egEvid_ipsi   = squeeze(nanmean(nanmean(egEvid(:,:,ipsiidx,midx),2),1));
          Xcontra       = repmat(xaxis(contraidx),[1 size(egEvid_contra,2)]);
          Xipsi         = repmat(xaxis(ipsiidx),[1 size(egEvid_ipsi,2)]);
          
          thisX          = Xcontra(:);
          thisY          = egEvid_contra(:);
          goodidx        = ~isnan(thisY);
          thisfit_contra = fit(abs(thisX(goodidx)),thisY(goodidx),'poly1');%,'robust','on');
          
          thisX          = Xipsi(:);
          thisY          = egEvid_ipsi(:);
          goodidx        = ~isnan(thisY);
          thisfit_ipsi   = fit(abs(thisX(goodidx)),thisY(goodidx),'poly1');%,'robust','on');

          ipsip1(iBoot)  = thisfit_ipsi.p1;
          contrap1(iBoot)= thisfit_contra.p1;
          ipsip2(iBoot)  = thisfit_ipsi.p2;
          contrap2(iBoot)= thisfit_contra.p2;
        end
        
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contraVSipsi.slope_diff_mean  = nanmean(contrap1 - ipsip1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contraVSipsi.slope_std        = nanstd(contrap1 - ipsip1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contraVSipsi.slope_isSig      = (nanmean(contrap1 - ipsip1) - cfg.nSDsig*nanstd(contrap1 - ipsip1) > 0) | (nanmean(contrap1 - ipsip1) + cfg.nSDsig*nanstd(contrap1 - ipsip1) < 0); 

        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.slope_mean  = nanmean(contrap1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.slope_std   = nanstd(contrap1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.slope_isSig = (nanmean(contrap1) - cfg.nSDsig*nanstd(contrap1) > 0) | (nanmean(contrap1) + cfg.nSDsig*nanstd(contrap1) < 0); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.offset_mean = nanmean(contrap2); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.offset_std  = nanstd(contrap2); 
        
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.slope_mean  = nanmean(ipsip1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.slope_std   = nanstd(ipsip1); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.slope_isSig = (nanmean(ipsip1) - cfg.nSDsig*nanstd(ipsip1) > 0) | (nanmean(ipsip1) + cfg.nSDsig*nanstd(ipsip1) < 0); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.offset_mean = nanmean(ipsip2); 
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.offset_std  = nanstd(ipsip2); 
        
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.xaxis = abs(xaxis(contraidx));
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.xaxis   = abs(xaxis(ipsiidx));
        
        egEvid_contra = squeeze(nanmean(nanmean(egEvid(:,:,contraidx,:),2),1));
        egEvid_ipsi   = squeeze(nanmean(nanmean(egEvid(:,:,ipsiidx,:),2),1));
        
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.dff_mean = mean(egEvid_contra,2);
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).contra.dff_sem  = std(egEvid_contra,0,2)./sqrt(nmice-1);
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.dff_mean   = mean(egEvid_ipsi,2);
        avgDffvsEvidence.linFits.(task).(trialType).mazeSection(iSection).ROI(iROI).ipsi.dff_sem    = std(egEvid_ipsi,0,2)./sqrt(nmice-1);
      end
    end
  end
end

end

%% ----------------------------------------------
%% plot fits
function [iPanel,figstats] = pnlFits(fig,iPanel,avgDffvsEvidence,cfg,egID,figstats)

% get relevant example
section    = cfg.egFitMzSect{egID};
iSection   = cellfun(@(x)(sum(x==section)==2),avgDffvsEvidence.linFits.mazeSections);
iROI       = strcmpi(avgDffvsEvidence.linFits.ROIlbl,cfg.egAreaFits{egID});
contra_dff = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).contra.dff_mean;
contra_sem = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).contra.dff_sem;
ipsi_dff   = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).ipsi.dff_mean;
ipsi_sem   = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).ipsi.dff_sem;

contra_a   = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).contra.slope_mean; 
contra_b   = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).contra.offset_mean; 
contraX    = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).contra.xaxis;
contraXfit = 0:15;
contraY    = contraXfit.*contra_a + contra_b;
% contra_aSD = avgDffvsEvidence.linFits.(cfg.egFitTask).(cfg.egFitTrialType).mazeSection(iSection).ROI(iROI).contra.slope_std;
% contra_bSD = avgDffvsEvidence.linFits.(cfg.egFitTask).(cfg.egFitTrialType).mazeSection(iSection).ROI(iROI).contra.offset_std;
ipsi_a     = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).ipsi.slope_mean; 
ipsi_b     = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).ipsi.offset_mean; 
ipsiX      = avgDffvsEvidence.linFits.(cfg.egFitTask{egID}).(cfg.egFitTrialType{egID}).mazeSection(iSection).ROI(iROI).ipsi.xaxis;
ipsiY      = contraXfit.*ipsi_a + ipsi_b;
% ipsi_aSD   = avgDffvsEvidence.linFits.(cfg.egFitTask).(cfg.egFitTrialType).mazeSection(iSection).ROI(iROI).ipsi.slope_std;
% ipsi_bSD   = avgDffvsEvidence.linFits.(cfg.egFitTask).(cfg.egFitTrialType).mazeSection(iSection).ROI(iROI).ipsi.offset_std;

%%
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

errorbar(contraX, contra_dff, contra_sem,'o', 'color', cfg.contraCl, 'markerfacecolor', cfg.contraCl);
errorbar(ipsiX, ipsi_dff, ipsi_sem,'o', 'color', cfg.ipsiCl, 'markerfacecolor', cfg.ipsiCl);
  
plot(contraXfit, contraY, '-', 'color', cfg.contraCl, 'linewidth', .75)
plot(contraXfit, ipsiY, '-', 'color', cfg.ipsiCl, 'linewidth', .75)

set(axs,'xtick',0:5:15); xlim([0 15])
xlabel('\Delta towers')
ylabel('\DeltaF/F (z-score)')
title(sprintf('%s, %s, %s\n%d - %d cm',cfg.egFitTask{egID},cfg.egFitTrialType{egID},cfg.egAreaFits{egID},section(1),section(2)))

yl = get(axs,'ylim');
text(1,yl(2)*.95,sprintf('Contra (slope = %1.2g)',contra_a),'color',cfg.contraCl)
text(1,yl(2)*.8,sprintf('Ipsi (slope = %1.2g)',ipsi_a),'color',cfg.ipsiCl)

figstats.panelA_right.xaxis_data          = contraX;
figstats.panelA_right.xaxis_fit           = contraXfit;
thislbl                                   = [cfg.egAreaFits{egID} '_' cfg.egFitTrialType{egID} 'Trials_' cfg.egFitTask{egID}];
figstats.panelA_right.(['mean_contra_' thislbl]) = contra_dff;
figstats.panelA_right.(['std_contra_' thislbl])  = contra_sem;
figstats.panelA_right.(['fit_contra_' thislbl])  = contraY;
figstats.panelA_right.(['mean_ipsi_' thislbl])   = ipsi_dff;
figstats.panelA_right.(['std_ipsi_' thislbl])    = ipsi_sem;
figstats.panelA_right.(['fit_ipsi_' thislbl])    = ipsiY;
figstats.panelA_right.xlabel             = '\Delta towers';
figstats.panelA_right.ylabel             = '\DeltaF/F (z-score)';


end

%% ----------------------------------------------
%% plot slope summary
function [iPanel,figstats] = pnlFitSlopeSummary(fig,iPanel,avgDffvsEvidence,cfg,task,trialtype,section,plotLegend,figstats)

%%
iSection    = cellfun(@(x)(sum(x==section)==2),avgDffvsEvidence.linFits.mazeSections);
lbls        = avgDffvsEvidence.linFits.ROIlbl;
nROI        = numel(lbls);
for iROI = 1:nROI
  contra_a(iROI)    = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).contra.slope_mean;
  contra_aSD(iROI)  = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).contra.slope_std;
  contra_aSig(iROI) = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).contra.slope_isSig;
  ipsi_a(iROI)      = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).ipsi.slope_mean;
  ipsi_aSD(iROI)    = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).ipsi.slope_std;
  ipsi_aSig(iROI)   = avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).ipsi.slope_isSig;
end

%%
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

xaxis = 1:3:nROI*3+1-3;

h(2) = bar(xaxis,contra_a,'edgecolor',cfg.contraCl,'facecolor','none','barwidth',.3);
h(4) = bar(xaxis+1,ipsi_a,'edgecolor',cfg.ipsiCl,'facecolor','none','barwidth',.3);

h(1) = bar(xaxis(contra_aSig),contra_a(contra_aSig),'edgecolor',cfg.contraCl,'facecolor',cfg.contraCl,'barwidth',.3);
h(3) = bar(xaxis(ipsi_aSig)+1,ipsi_a(ipsi_aSig),'edgecolor',cfg.ipsiCl,'facecolor',cfg.ipsiCl,'barwidth',.3);

errorbar(xaxis,contra_a,contra_aSD,'.','markersize',0.1,'color',cfg.contraCl);
errorbar(xaxis+1,ipsi_a,ipsi_aSD,'.','markersize',0.1,'color',cfg.ipsiCl);

% contra vs ipsi significance
for iROI = 1:nROI
  if avgDffvsEvidence.linFits.(task).(trialtype).mazeSection(iSection).ROI(iROI).contraVSipsi.slope_isSig
    text(xaxis(iROI)+.5,-.25,'*','fontsize',12,'color','k')
  end
end

xlim([0 xaxis(end)+2])
set(axs,'xtick',xaxis+.5,'xticklabel',lbls)
ylabel('Slope (\DeltaF/F/\Delta towers)(z-score)')
if plotLegend; legend(h,{'contra (sig.)','contra (n.s.)','ipsi (sig.)','ipsi (n.s.)'}); end
set(axs,'view',[90 -90])
title(sprintf('%s, %s\n%d - %d cm',task,trialtype,section(1),section(2)))


switch task
  case 'accumul'
    thislbl                            = [trialtype 'Trials_' task];
    figstats.panelB.(['contra_slope_avg_' thislbl])           = contra_a;
    figstats.panelB.(['contra_slope_std_' thislbl])           = contra_aSD;
    figstats.panelB.(['contra_slope_isSignificant_' thislbl]) = contra_aSig;
    figstats.panelB.(['ipsi_slope_avg_' thislbl])             = ipsi_a;
    figstats.panelB.(['ipsi_slope_std_' thislbl])             = ipsi_aSD;
    figstats.panelB.(['ipsi_slope_isSignificant_' thislbl])   = ipsi_aSig;
    figstats.panelB.ylabel             = 'Slope (\DeltaF/F/\Delta towers)(z-score)';
    figstats.panelB.ROIlabels          = lbls;
  otherwise
    thislbl                            = [trialtype 'Trials_' task];
    figstats.panelC.(['contra_slope_avg_' thislbl])           = contra_a;
    figstats.panelC.(['contra_slope_std_' thislbl])           = contra_aSD;
    figstats.panelC.(['contra_slope_isSignificant_' thislbl]) = contra_aSig;
    figstats.panelC.(['ipsi_slope_avg_' thislbl])             = ipsi_a;
    figstats.panelC.(['ipsi_slope_std_' thislbl])             = ipsi_aSD;
    figstats.panelC.(['ipsi_slope_isSignificant_' thislbl])   = ipsi_aSig;
    figstats.panelC.ylabel             = 'Slope (\DeltaF/F/\Delta towers)(z-score)';
    figstats.panelC.ROIlabels          = lbls;
end
end

% %% ----------------------------------------------
% %% plot avg Dff per evidence
% function [iPanel,evidTunSumm] = pnlGPtuning(fig, iPanel, evidTunSumm, cfg, task)
% 
% % make it into contra/ipsi, then do contra - ipsi
% yhatmat    = evidTunSumm.(cfg.GPevType).(cfg.GPepoch).yhat;
% xaxis      = evidTunSumm.(cfg.GPevType).(cfg.GPepoch).xaxis;
% yhatmatL   = yhatmat(:,1:2:end,:);
% yhatmatR   = flipud(yhatmat(:,2:2:end,:));
% yhatmat    = (yhatmatL+yhatmatR)./2;
% 
% nEv           = numel(xaxis);
% yhatmatContra = yhatmat((nEv-1)/2+2:nEv,:,:);
% yhatmatIpsi   = yhatmat(1:(nEv-1)/2,:,:);
% yhatmat       = yhatmatContra - yhatmatIpsi;
% xaxis         = 1:(nEv-1)/2;
% evidTunSumm.yhat_contraMinusIpsi  = yhatmat;
% evidTunSumm.xaxis_contraMinusIpsi = xaxis;
% 
% % ANOVA to compare ROIs
% [nev,nROI,nmice] = size(yhatmat);
% evidMat          = repmat(xaxis',[1 nROI nmice]);
% ROIMat           = repmat(1:nROI,[nev 1 nmice]);
% mouseMat         = repmat(reshape(1:nmice,[1 1 nmice]),[nev nROI 1]);
% [evidTunSumm.stats.ANOVA_p,evidTunSumm.stats.ANOVA_table,evidTunSumm.stats.ANOVA_stats] ...
%                  = anovan(yhatmat(:),{evidMat(:),ROIMat(:),mouseMat(:)},'display','off','varnames',{'evidence','ROI','mouse'});
% 
% % plot
% nROI     = size(yhatmat,2);
% lbls     = unique(cellfun(@(x)(x(1:end-2)),evidTunSumm.ROIlbl,'UniformOutput',false),'stable');
% cl       = getDefaultROIcl(lbls);
% 
% iPanel                  = iPanel + 1;
% axs                     = fig.panel(iPanel);
% hold(axs,'on')
% 
% for iROI = 1:nROI 
%   plot(xaxis, nanmean(yhatmat(:,iROI,:),3),'-', 'LineWidth', .75, 'color', cl(iROI,:));
% end
% 
% xlabel('\Delta towers')
% ylabel(sprintf('%sF/F (z-score)\n(contra - ipsi towers)','\Delta'))
% title(task)
% 
% xlim([1 10]); 
% if strcmp(task,'Accum. towers')
%   ylim([-.15 .25])
% else
%   ylim([-.3 .5])
% end
% 
% % text(2,.2,sprintf('p(ROI) = %1.2g',evidTunSumm.stats.ANOVA_p(2)),'fontsize',9,'color','k')
% end

% %% ----------------------------------------------
% %% plot avg Dff per evidence
% function [iPanel,evidTunSumm] = pnlTaskComp(fig, iPanel, evidTunSumm)
% 
% % max - min evidence for both tasks
% yhat = evidTunSumm.yhat_contraMinusIpsi;
% 
% % stats
% 
% % plot
% end