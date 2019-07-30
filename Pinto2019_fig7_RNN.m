function RNNdata = Pinto2019_fig7_RNN(analysisFilePath,summaryFile)

% RNNdata = Pinto2019_fig7_RNN(analysisFilePath,summaryFile)
% plots RNN results, except network schematics

%% ------------------------------------------------------------------------
%% Instructions for running analysis 
% 1) To train network and obtain correlations: RNN.m (repo: behaviorAnalysis)
% 2) Inactivations: inactivate.m and inactivate_module.m
% 3) example use: RNNdemo.m
% 4) summary data: summarizRNNbyPconn.m, summarizRNNinactivationByPconn.m
% 5) summary data used to plot figures are saved on braininit bucket, raw
%    data files from simulations can be found in //bucket/brody/briandd/RNN/data/
% RNN code, except summary data, was written by Brian DePasquale & Kanaka Rajan

%%
cfg.clustCl       = [60 179 113; 10 35 140]./255; % module colors
cfg.histBins      = -1:.1:1; % for cc histogram
cfg.withinClustCl = [230 186 0]./255; 
cfg.xClustCl      = [.6 .6 .6];
cfg.towersCl      = [0 0 0]; 
cfg.ctrlCl        = [.6 .6 .6]; 
cfg.visEgIdx      = [110 113 490]; % unit IDs for rate examples
cfg.frontEgIdx    = [604 690 700]; % unit IDs for rate examples
cfg.PconnLevel    = 0.05; % level of Pcon used for main RNN
cfg.plotwhat      = 'rate'; % which inactivation measure
cfg.sameG         = false; % load files with same g in on- and off-diagonals

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
if cfg.sameG
  ext = '_sameG.mat';
else
  ext = '.mat';
end
ratefn   = sprintf('%sdata_and_network_sparse_%g%s',analysisFilePath,cfg.PconnLevel,ext);
inactfn  = [analysisFilePath 'RNNinactivationSummaryByPconn' ext];

rates    = load(ratefn,'post_post','ant_post','Task_data','N','R2data','Rdata_untrained','T','dt','N');
inact    = load(inactfn);

%% correlations between rates of the RNN units
visIdx      = rates.post_post;
frontIdx    = rates.ant_post;

pulsesR     = rates.R2data([visIdx frontIdx],:,rates.Task_data == 1);
detectR     = rates.R2data([visIdx frontIdx],:,rates.Task_data == 2);

nPulseTrial = size(pulsesR,3);
EA_corr     = zeros(rates.N,rates.N,nPulseTrial);
for iTrial = 1:nPulseTrial
    EA_corr(:,:,iTrial) = corr(pulsesR(:,:,iTrial)');
end

nDetectTrial = size(detectR,3);
detect_corr  = zeros(rates.N,rates.N,nDetectTrial);
for iTrial = 1:nDetectTrial
    detect_corr(:,:,iTrial) = corr(detectR(:,:,iTrial)');
end

EA_corr     = mean(EA_corr,3);
detect_corr = mean(detect_corr,3);
corrDiff    = EA_corr - detect_corr;
corrMagDiff = abs(EA_corr) - abs(detect_corr);

front       = corrMagDiff(frontIdx,frontIdx);
vis         = corrMagDiff(visIdx,visIdx);
xClust      = corrMagDiff(visIdx,frontIdx);

xClust      = nanmean(xClust,2);
withinClust = mean([nanmean(front,2) nanmean(vis,2)],2);

%% compile data and stats 
RNNdata.rates.untrainedEg   = rates.Rdata_untrained;
RNNdata.rates.accumulTaskEg = pulsesR(:,:,end);
RNNdata.rates.ctrlTaskEg    = detectR(:,:,end);
RNNdata.rates.visIdx        = visIdx;
RNNdata.rates.frontIdx      = frontIdx;
RNNdata.rates.absCCdiff     = corrMagDiff;
RNNdata.rates.CCdiff        = corrDiff;

% corr stats
cc                              = abs(EA_corr);
cc(logical(eye(size(cc))))      = nan;
RNNdata.rates.avgCCunit_accumul = nanmean(cc,2);
cc                              = abs(detect_corr);
cc(logical(eye(size(cc))))      = nan;
RNNdata.rates.avgCCunit_ctrl    = nanmean(cc,2);
RNNdata.rates.avgCCunit_pval    = signrank(RNNdata.rates.avgCCunit_accumul,RNNdata.rates.avgCCunit_ctrl);
RNNdata.rates.p_ccDelta_withinVSacrossModule = signrank(withinClust,xClust);

%% inactivation stats

% random
idx                             = inact.PconnVals == cfg.PconnLevel;
RNNdata.inact_rand.EA_perf      = nanmean(inact.EA_perf{idx});
RNNdata.inact_rand.EA_perf_sem  = nanstd(inact.EA_perf{idx})./sqrt(inact.nruns-1);
RNNdata.inact_rand.EA_delta     = inact.EA_fitParams{idx}(:,2);
RNNdata.inact_rand.EA_rate      = inact.EA_fitParams{idx}(:,3);
RNNdata.inact_rand.EA_tolerance = inact.EA_tolerance_fit(:,idx);

RNNdata.inact_rand.detect_perf      = nanmean(inact.detect_perf{idx});
RNNdata.inact_rand.detect_perf_sem  = nanstd(inact.detect_perf{idx})./sqrt(inact.nruns-1);
RNNdata.inact_rand.detect_delta     = inact.detect_fitParams{idx}(:,2);
RNNdata.inact_rand.detect_rate      = inact.detect_fitParams{idx}(:,3);
RNNdata.inact_rand.detect_tolerance = inact.detect_tolerance_fit(:,idx);

% random, task comparison
RNNdata.inact_rand.p_delta      = signrank(RNNdata.inact_rand.EA_delta,RNNdata.inact_rand.detect_delta);
RNNdata.inact_rand.p_rate       = signrank(RNNdata.inact_rand.EA_rate,RNNdata.inact_rand.detect_rate);
RNNdata.inact_rand.p_tolerance  = signrank(RNNdata.inact_rand.EA_tolerance,RNNdata.inact_rand.detect_tolerance);

% module
RNNdata.inact_vis.EA_perf      = mean(inact.EA_perf_post);
RNNdata.inact_vis.EA_perf_sem  = std(inact.EA_perf_post)./sqrt(inact.nruns-1);
RNNdata.inact_vis.EA_delta     = inact.EA_fitParams_post(:,2);
RNNdata.inact_vis.EA_rate      = inact.EA_fitParams_post(:,3);
RNNdata.inact_vis.EA_tolerance = inact.EA_tolerance_fit_post;

RNNdata.inact_vis.detect_perf      = mean(inact.detect_perf_post);
RNNdata.inact_vis.detect_perf_sem  = std(inact.detect_perf_post)./sqrt(inact.nruns-1);
RNNdata.inact_vis.detect_delta     = inact.detect_fitParams_post(:,2);
RNNdata.inact_vis.detect_rate      = inact.detect_fitParams_post(:,3);
RNNdata.inact_vis.detect_tolerance = inact.detect_tolerance_fit_post;

RNNdata.inact_front.EA_perf      = mean(inact.EA_perf_ant);
RNNdata.inact_front.EA_perf_sem  = std(inact.EA_perf_ant)./sqrt(inact.nruns-1);
RNNdata.inact_front.EA_delta     = inact.EA_fitParams_ant(:,2);
RNNdata.inact_front.EA_rate      = inact.EA_fitParams_ant(:,3);
RNNdata.inact_front.EA_tolerance = inact.EA_tolerance_fit_ant;

RNNdata.inact_front.detect_perf      = mean(inact.detect_perf_ant);
RNNdata.inact_front.detect_perf_sem  = std(inact.detect_perf_ant)./sqrt(inact.nruns-1);
RNNdata.inact_front.detect_delta     = inact.detect_fitParams_ant(:,2);
RNNdata.inact_front.detect_rate      = inact.detect_fitParams_ant(:,3);
RNNdata.inact_front.detect_tolerance = inact.detect_tolerance_fit_ant;

% ANOVA for inactivation level and task
perfEA  = inact.EA_perf{idx}(:,1:end-1);
perfD   = inact.detect_perf{idx}(:,1:end-1);
fracMat = ones(inact.nruns,inact.nfractions-1)*diag(inact.percents(1:end-1)); 
perfvec = [perfEA(:); perfD(:)];
fracVec = repmat(fracMat(:),[2 1]);
taskvec = [ones(numel(perfEA),1) 2.*ones(numel(perfEA),1)];
isgood  = ~isnan(perfvec);

RNNdata.randInactANOVA = anovan(perfvec(isgood),{taskvec(isgood),fracVec(isgood)}, ...
                                'display','off','varnames',{'module','task'});

% ANOVA for module and task
deltavec = [RNNdata.inact_vis.EA_delta; RNNdata.inact_vis.detect_delta; RNNdata.inact_front.EA_delta; RNNdata.inact_front.detect_delta];
ratevec  = [RNNdata.inact_vis.EA_rate; RNNdata.inact_vis.detect_rate; RNNdata.inact_front.EA_rate; RNNdata.inact_front.detect_rate];
tolvec   = [RNNdata.inact_vis.EA_tolerance; RNNdata.inact_vis.detect_tolerance; RNNdata.inact_front.EA_tolerance; RNNdata.inact_front.detect_tolerance];
modvec   = [ones(inact.nruns*2,1);2.*ones(inact.nruns*2,1)];
taskvec  = repmat([ones(inact.nruns,1);2.*ones(inact.nruns,1)],[2 1]);


[RNNdata.moduleInactANOVA_delta,~,RNNdata.moduleInactANOVAstats_delta] = ...
  anovan(deltavec,{modvec;taskvec},'display','off','varnames',{'module','task'});
[RNNdata.moduleInactANOVA_rate,~,RNNdata.moduleInactANOVAstats_rate] = ...
  anovan(ratevec,{modvec;taskvec},'display','off','varnames',{'module','task'});
[RNNdata.moduleInactANOVA_tolerance,~,RNNdata.moduleInactANOVAstats_tolerance] = ...
  anovan(tolvec,{modvec;taskvec},'display','off','varnames',{'module','task'});

perfEA_ant  = inact.EA_perf_ant(:,1:end-1);
perfD_ant   = inact.detect_perf_ant(:,1:end-1);
perfEA_post = inact.EA_perf_post(:,1:end-1);
perfD_post  = inact.detect_perf_post(:,1:end-1);
perfvec     = [perfEA_ant(:); perfD_ant(:);perfEA_post(:); perfD_post(:)];
fracVec     = repmat(fracMat(:),[4 1]);
taskvec     = repmat([ones(numel(perfEA_ant),1); 2.*ones(numel(perfEA_ant),1)],[2 1]);
modvec      = [ones(numel(perfEA_ant)*2,1); 2.*ones(numel(perfEA_ant)*2,1)];

[RNNdata.moduleInactANOVA_perf, ~ ,RNNdata.moduleInactANOVAstats_perf] = ...
  anovan(perfvec,{modvec;taskvec;fracVec},'display','off','varnames',{'module','task','percent'});

%% FIGURE
layout       = [   1  2 13 13 14 14 ...
                ;  5  6 13 13 14 14 ...
                ;  9 10 15 15 16 16 ...
                ;  3  4 15 15 16 16 ...
                ;  7  8 17 17 17 17 ...
                ; 11 12 18 18 18 18 ...
               ];
fig          = PaneledFigure(layout, 'tiny');
iPanel       = 0;


%% examples
xaxis          = 0:rates.dt:rates.T;
fig7stats.panelB.xaxis = xaxis;

for iEg = 1:numel(cfg.visEgIdx)
  
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  plot([0 300],[0 0],'--','color',[.7 .7 .7],'linewidth',.25)
  plot(xaxis,RNNdata.rates.ctrlTaskEg(cfg.frontEgIdx(iEg),:),'-','color',cfg.clustCl(2,:),'linewidth',1);
  box off; ylim([-1 1]); xlim([0 300]); set(gca,'xtick',[],'xcolor','w')
  if iEg == 1; title('guided'); end
  
  fig7stats.panelB.(['frontalModule_visGuided_eg' num2str(iEg)]) = RNNdata.rates.ctrlTaskEg(cfg.frontEgIdx(iEg),:);
  
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  plot([0 300],[0 0],'--','color',[.7 .7 .7],'linewidth',.25)
  plot(xaxis,RNNdata.rates.accumulTaskEg(cfg.frontEgIdx(iEg),:),'-','color',cfg.clustCl(2,:),'linewidth',1);
  box off; ylim([-1 1]); xlim([0 300]); set(gca,'xtick',[],'xcolor','w')
  if iEg == 1; title('accumul'); end
  
  fig7stats.panelB.(['frontalModule_accumTowers_eg' num2str(iEg)]) = RNNdata.rates.accumulTaskEg(cfg.frontEgIdx(iEg),:);
  
  
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  plot([0 300],[0 0],'--','color',[.7 .7 .7],'linewidth',.25)
  plot(xaxis,RNNdata.rates.ctrlTaskEg(cfg.visEgIdx(iEg),:),'-','color',cfg.clustCl(1,:),'linewidth',1);
  box off; ylim([-1 1]); xlim([0 300]); set(gca,'xtick',[],'xcolor','w')
  
  fig7stats.panelB.(['posteriorModule_visGuided_eg' num2str(iEg)]) = RNNdata.rates.ctrlTaskEg(cfg.visEgIdx(iEg),:);
  
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  plot([0 300],[0 0],'--','color',[.7 .7 .7],'linewidth',.25)
  plot(xaxis,RNNdata.rates.accumulTaskEg(cfg.visEgIdx(iEg),:),'-','color',cfg.clustCl(1,:),'linewidth',1);
  box off; ylim([-1 1]); xlim([0 300]); set(gca,'xtick',[],'xcolor','w')
  
  fig7stats.panelB.(['posteriorModule_visGuided_eg' num2str(iEg)]) = RNNdata.rates.accumulTaskEg(cfg.visEgIdx(iEg),:);
end

fig7stats.panelB.xlabel = 'Time (s)';
fig7stats.panelB.ylabel = 'Rate (a.u.)';

%% corr matrices
visIdx      = [.5 500.5];
frontIdx    = [500.5 1000.5];

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

imagesc(abs(EA_corr),[0 1]); colormap gray; axis ij; axis tight
plot([visIdx(1) visIdx(end) visIdx(end) visIdx(1)   visIdx(1)], ...
     [visIdx(1) visIdx(1)   visIdx(end) visIdx(end) visIdx(1)], ...
     '-','color',cfg.clustCl(1,:),'linewidth',2  );
plot([frontIdx(1) frontIdx(end) frontIdx(end) frontIdx(1)   frontIdx(1)], ...
     [frontIdx(1) frontIdx(1)   frontIdx(end) frontIdx(end) frontIdx(1)], ...
     '-','color',cfg.clustCl(2,:),'linewidth',2 );
set(axs,'xtick',[],'ytick',[])
ylabel ('RNN units'); xlabel('RNN units'); title('Accumulation')

fig7stats.panelC.absCC_accumulTowers = abs(EA_corr);

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')
imagesc(abs(detect_corr),[0 1]); colormap gray; axis off; axis ij; axis tight
plot([visIdx(1) visIdx(end) visIdx(end) visIdx(1)   visIdx(1)], ...
     [visIdx(1) visIdx(1)   visIdx(end) visIdx(end) visIdx(1)], ...
     '-','color',cfg.clustCl(1,:),'linewidth',2  );
plot([frontIdx(1) frontIdx(end) frontIdx(end) frontIdx(1)   frontIdx(1)], ...
     [frontIdx(1) frontIdx(1)   frontIdx(end) frontIdx(end) frontIdx(1)], ...
     '-','color',cfg.clustCl(2,:),'linewidth',2  );
set(axs,'xtick',[],'ytick',[])
title('Control')

fig7stats.panelC.absCC_visGuided = abs(detect_corr);
fig7stats.panelC.dataLabel       = 'abs(r)';

cbar              = smallcolorbar(axs,'southoutside');
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = 'abs(r)';

%% histogram of task-induced corr changes
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot([0 0],[0 .5],'--','color',[.7 .7 .7])
xaxis                   = toBinCenters(cfg.histBins);
countsWithin            = histcounts(withinClust,cfg.histBins)./numel(withinClust);
countsAcross            = histcounts(xClust,cfg.histBins)./numel(xClust);
hHisto                  = gobjects(0);
hHisto(end+1)           = bar ( axs, xaxis, countsWithin, .96              ...
                              , 'EdgeColor'     , 'none'                   ...
                              , 'FaceColor'     , cfg.withinClustCl        ...
                              );
                            
                            
hHisto(end+1)           = bar ( axs, xaxis, countsAcross, .96              ...
                              , 'EdgeColor'     , 'none'                   ...
                              , 'FaceColor'     , cfg.xClustCl             ...
                              );

% Simulate transparency (actual transparency has PDF output issues)
colors                  = [cfg.withinClustCl; cfg.xClustCl];
                          bar ( axs, xaxis, min([countsWithin; countsAcross],[],1), .96     ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , mean(colors,1)                            ...
                              );
                            
legend(axs,hHisto([1 2]),{'within module','across module'})
xlabel(['\Delta abs(corr)' sprintf('\n(accumul - ctrl)')])
ylabel('Prop. units')

text(.3,.25,sprintf('P = %1.2g',RNNdata.rates.p_ccDelta_withinVSacrossModule))
xlim([-1 1]); ylim([0 .5])

fig7stats.panelD.deltaCC_withinModule = countsWithin;
fig7stats.panelD.deltaCC_acrossModule = countsAcross;
fig7stats.panelD.xaxis                = xaxis;
fig7stats.panelD.xlabel               = ['\Delta abs(corr)' sprintf('\n(accumul.-towers - vis.-guided)')];
fig7stats.panelD.ylabel               = 'Prop. units';

%% task performance (inactivation)
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot([0 inact.percents(end)],[50 50],'--','color',[.7 .7 .7])
h(1) = errorbar(inact.percents,RNNdata.inact_rand.detect_perf,RNNdata.inact_rand.detect_perf_sem,'.','color',cfg.ctrlCl,'markersize',1);
h(2) = errorbar(inact.percents,RNNdata.inact_rand.EA_perf,RNNdata.inact_rand.EA_perf_sem,'.','color',cfg.towersCl,'markersize',1);

fig7stats.panelE.inactivation_accumTowers_mean = RNNdata.inact_rand.EA_perf;
fig7stats.panelE.inactivation_accumTowers_sem  = RNNdata.inact_rand.EA_perf_sem;
fig7stats.panelE.inactivation_visGuided_mean   = RNNdata.inact_rand.detect_perf;
fig7stats.panelE.inactivation_visGuided_sem    = RNNdata.inact_rand.detect_perf_sem;
fig7stats.panelE.xaxis                         = inact.percents;

thisf = fit(inact.percents,RNNdata.inact_rand.detect_perf',inact.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inact.percents(end),fiteval(thisf,(0:.01:inact.percents(end))'),'-','color',cfg.ctrlCl,'linewidth',.75)

fig7stats.panelE.xaxis_exponentialFit            = 0:.01:inact.percents(end);
fig7stats.panelE.yaxis_exponentialFit_visGuided  = fiteval(thisf,(0:.01:inact.percents(end))');

thisf = fit(inact.percents,RNNdata.inact_rand.EA_perf',inact.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inact.percents(end),fiteval(thisf,(0:.01:inact.percents(end))'),'-','color',cfg.towersCl,'linewidth',.75)

fig7stats.panelE.yaxis_exponentialFit_accumTowers = fiteval(thisf,(0:.01:inact.percents(end))');

legend(h,{'guided','accum.'});
xlabel('Inactivated units (%)')
ylabel('Perf. (% correct)')
% xlim([0 inact.percents(end)])
xlim([0 10])
ylim([45 100])

if contains(cfg.plotwhat,'tolerance')
  vec1 = RNNdata.(['inact_rand.detect_' cfg.plotwhat]);
  vec2 = RNNdata.(['inact_rand.EA_' cfg.plotwhat]);
  plot(mean(vec1),'v','color',cfg.ctrlCl)
  plot(mean(vec2),'v','color',cfg.towersCl)
end

fig7stats.panelE.xlabel               = 'Inactivated units (%)';
fig7stats.panelE.ylabel               = 'Perf. (% correct)';


%% Save data and figure to disk
save(summaryFile,'RNNdata','fig7stats','-append')

fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])