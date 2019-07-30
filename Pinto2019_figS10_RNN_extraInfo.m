function Pinto2019_figS10_RNN_extraInfo(analysisFilePath,summaryFile)

% RNNdata = Pinto2019_figS10_RNN_extraInfo(analysisFilePath,summaryFile)
% plots supplemental details on RNN results

%% ------------------------------------------------------------------------
%% Instructions for running analysis 
% 1) To train network and obtain correlations: RNN.m ( behaviorAnalysis)
% 2) Inactivations: inactivate.m and inactivate_module.m
% 3) example use: RNNdemo.m
% 4) summary data: summarizRNNbyPconn.m, summarizRNNinactivationByPconn.m
% RNN code, except summary data, was written by Brian DePasquale & Kanaka Rajan

%%
cfg.clustCl              = [60 179 113; 10 35 140]./255; % module colors
cfg.withinClustCl        = [230 186 0]./255;
cfg.xClustCl             = [.6 .6 .6];
cfg.towersCl             = [0 0 0]; 
cfg.ctrlCl               = [.6 .6 .6]; 
cfg.decoderEgTrialAccum  = [4 1]; % decoder egs, will only use first
cfg.decoderEgTrialGuide  = [8 9]; % decoder egs, will only use first
cfg.PconnVals            = [.05 .1 .25 .5 1]; % tested teacher Pconn
cfg.normByg              = false; % to normalize off-diag and diag blocks by their gs
cfg.matTh                = true; % to threshold weight mat
cfg.plotwhat             = 'rate'; % inactivation quantification
cfg.sameG                = false; % use files with same g on on- and off-diagonals
cfg.PconnLevel           = 0.05; % Pconn level used in main figure
cfg.weightClim           = [-.2 .2]; % limit for weight matrices
cfg.nIterBoot            = 100; % bootstrapping iterations for exponential fit
cfg.model                = 'a+b*exp(c*x)'; % fit perf decay
cfg.r2th                 = 0.7; % minimal fit r2 
cfg.globalSigStrength    = .75; % how much of a global signal to add for correlations
cfg.rbins                = -1:.05:1; % for correlation histograms

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});

%% load data
decoderfn          = sprintf('%sdata_and_network_sparse_%g.mat',analysisFilePath,cfg.PconnLevel);
decoder            = load(decoderfn,'zsdata','outputdata','pulsedata');
corrByConn         = load([analysisFilePath 'RNNsummaryByPconn.mat']);
inactByConn        = load([analysisFilePath 'RNNinactivationSummaryByPconn.mat']);

%% correlations between rates of the RNN units
xClustCorr_accum_mean = zeros(1,numel(cfg.PconnVals));
xClustCorr_accum_std  = zeros(1,numel(cfg.PconnVals));
xClustCorr_guide_mean = zeros(1,numel(cfg.PconnVals));
xClustCorr_guide_std  = zeros(1,numel(cfg.PconnVals));

for iConn = 1:numel(cfg.PconnVals)
  thiscorr = abs(corrByConn.EA_corr{iConn}(corrByConn.visIdx,corrByConn.frontIdx));
  xClustCorr_accum_mean(iConn) = mean(mean(thiscorr,2));
  xClustCorr_accum_std(iConn)  = std(mean(thiscorr,2));
  
  thiscorr = abs(corrByConn.detect_corr{iConn}(corrByConn.visIdx,corrByConn.frontIdx));
  xClustCorr_guide_mean(iConn) = mean(mean(thiscorr,2));
  xClustCorr_guide_std(iConn)  = std(mean(thiscorr,2));
end

%% FIGURE
layout       = [   1  3  5  7  9 ...
                ;  2  4  6  8 10 ...
                ; 11 11 13 13 15 ...
                ; 12 12 14 14 16 ...
                ; 17 18 19 20 21 ...
               ];
fig          = PaneledFigure(layout, 'tiny');
iPanel       = 0;

%% J for teacher and learner
visIdx   = corrByConn.visIdx+.5;
frontIdx = corrByConn.frontIdx+.5;
for iConn = 1:numel(cfg.PconnVals)
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  if cfg.normByg
    mat   = corrByConn.teacherJ{iConn};
    mat(1:500,1:500)       = mat(1:500,1:500)./corrByConn.gvals(iConn);
    mat(501:1000,501:1000) = mat(501:1000,501:1000)./corrByConn.gvals(iConn);
    mat(1:500,501:1000)    = mat(1:500,501:1000)./corrByConn.goffvals(iConn);
    mat(501:1000,1:500)    = mat(501:1000,1:500)./corrByConn.goffvals(iConn);
    clims = [-.25 .25];
  elseif cfg.matTh
    mat   = corrByConn.teacherJ{iConn};
    mat(mat<0) = -1;
    mat(mat>0) = 1;
    clims = [-1 1];
  else
    mat   = corrByConn.teacherJ{iConn};
    clims = [-.5 .5];
  end
  imagesc(mat,clims); colormap red2blue; axis ij; axis tight
  plot([visIdx(1) visIdx(end) visIdx(end) visIdx(1)   visIdx(1)], ...
       [visIdx(1) visIdx(1)   visIdx(end) visIdx(end) visIdx(1)], ...
       '-','color',cfg.clustCl(1,:),'linewidth',2  );
  plot([frontIdx(1) frontIdx(end) frontIdx(end) frontIdx(1)   frontIdx(1)], ...
       [frontIdx(1) frontIdx(1)   frontIdx(end) frontIdx(end) frontIdx(1)], ...
       '-','color',cfg.clustCl(2,:),'linewidth',2 );
  set(axs,'xtick',[],'ytick',[])
  if iConn == 1; ylabel ('Teacher units'); xlabel('Teacher units'); end
  title(['P_{conn} = ' num2str(cfg.PconnVals(iConn))],'fontsize',12)
  box on
  if iConn == 5
    cbar              = smallcolorbar(axs,'southoutside');
    pos               = cbar.Position;
    cbar.Units        = get(axs,'units');
    cbar.Position     = pos;
    cbar.Label.String = 'weight sign';
  end
  
  iPanel      = iPanel + 1;
  axs         = fig.panel(iPanel);
  hold(axs, 'on')
  if cfg.normByg
    mat = corrByConn.learnerJ{iConn}./corrByConn.gvals(iConn);
  else
    mat   = corrByConn.learnerJ{iConn};
  end
  imagesc(mat,cfg.weightClim); colormap red2blue; axis ij; axis tight
  plot([visIdx(1) visIdx(end) visIdx(end) visIdx(1)   visIdx(1)], ...
       [visIdx(1) visIdx(1)   visIdx(end) visIdx(end) visIdx(1)], ...
       '-','color',cfg.clustCl(1,:),'linewidth',2  );
  plot([frontIdx(1) frontIdx(end) frontIdx(end) frontIdx(1)   frontIdx(1)], ...
       [frontIdx(1) frontIdx(1)   frontIdx(end) frontIdx(end) frontIdx(1)], ...
       '-','color',cfg.clustCl(2,:),'linewidth',2 );
  set(axs,'xtick',[],'ytick',[])
  if iConn == 1; ylabel ('Learner units'); xlabel('Learner units'); end
  if iConn == 5
    cbar              = smallcolorbar(axs,'southoutside');
    pos               = cbar.Position;
    cbar.Units        = get(axs,'units');
    cbar.Position     = pos;
    cbar.Label.String = 'weight (a.u.)';
  end
  box on
end

%% example decoding, accumul
pulsesR = decoder.pulsedata(1,:,cfg.decoderEgTrialAccum(1));
pulsesL = decoder.pulsedata(2,:,cfg.decoderEgTrialAccum(1));
decoded = decoder.zsdata(2,:,cfg.decoderEgTrialAccum(1));
actual  = decoder.outputdata(2,:,cfg.decoderEgTrialAccum(1));
time    = linspace(0,3,numel(decoded));

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot(time,.3+pulsesR,'-','color',widefieldParams.myblue)
plot(time,pulsesL,'-','color',widefieldParams.myred)
set(axs,'xtick',[],'ytick',[])
title('Accumul task')

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot(time,actual,'-','color',[.5 .5 .5])
plot(time,decoded,'-','color','k')
xlabel('time (s)')
set(axs,'xtick',0:1:3)
ylabel(sprintf('Proj. onto\n accumul. coding axis'))
legend({'actual','decoded'},'location','best')
ylim([-1.3 1.3])
set(axs,'ytick',-1:1:1)

%% example decoding, guided
pulsesR = decoder.pulsedata(1,:,cfg.decoderEgTrialGuide(1));
pulsesL = decoder.pulsedata(2,:,cfg.decoderEgTrialGuide(1));
decoded = decoder.zsdata(2,:,cfg.decoderEgTrialGuide(1));
actual  = decoder.outputdata(2,:,cfg.decoderEgTrialGuide(1));
time    = linspace(0,3,numel(decoded));

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot(time,.3+pulsesR,'-','color',widefieldParams.myblue)
plot(time,pulsesL,'-','color',widefieldParams.myred)
if sum(pulsesR) > sum(pulsesL)
  plot([time(1) time(end)],[.3+max(pulsesR) .3+max(pulsesR)],'-','color',widefieldParams.myblueSh)
else
  plot([time(1) time(end)],[max(pulsesL) max(pulsesL)],'-','color',widefieldParams.myredSh)
end
set(axs,'xtick',[],'ytick',[])
title('Guided task')

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot(time,actual,'-','color',[.5 .5 .5])
plot(time,decoded,'-','color','k')
xlabel('time (s)')
set(axs,'xtick',0:1:3)
ylabel(sprintf('Proj. onto\n accumul. coding axis'))
legend({'actual','decoded'},'location','best')
ylim([-1.3 1.3])
set(axs,'ytick',-1:1:1)

%% x-clust corr by Pconn
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

errorbar(cfg.PconnVals,xClustCorr_accum_mean,xClustCorr_accum_std, ...
         'color',cfg.towersCl,'linewidth',.75)
errorbar(cfg.PconnVals,xClustCorr_guide_mean,xClustCorr_guide_std, ...
         'color',cfg.ctrlCl,'linewidth',.75)
legend({'accum','guided'},'location','best')
set(axs,'xscale','log','xtick',[.05 .5 1])
xlabel('Teacher P_{conn}')
ylabel('Inter-module corr (|r|)')
xlim([.025 1.5])
ylim([.2 .8])

%% inactivation by Pconn

% for each Pconn fit exponential to the mean and bootstrap for error
EA_rate     = nan(cfg.nIterBoot,numel(cfg.PconnVals));
detect_rate = EA_rate;
for iConn = 1:numel(cfg.PconnVals)
  EA_perf     = inactByConn.EA_perf{iConn};
  detect_perf = inactByConn.detect_perf{iConn};
  nruns       = size(EA_perf,1);
  for iBoot = 1:cfg.nIterBoot
    idx    = randsample(nruns,nruns,true);
    EAmean = mean(EA_perf(idx,:))';
    DTmean = mean(detect_perf(idx,:))';
    [thisf,stats]  = fit(inactByConn.percents,EAmean,cfg.model,'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
    if stats.rsquare > cfg.r2th; EA_rate(iBoot,iConn) = thisf.c; end
    [thisf,stats]  = fit(inactByConn.percents,DTmean,cfg.model,'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
    if stats.rsquare > cfg.r2th; detect_rate(iBoot,iConn) = thisf.c; end
  end
end
%%
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

RNNdataByPconn.EA_inact_mean     = nanmean(-EA_rate);
RNNdataByPconn.EA_inact_std      = nanstd(-EA_rate);
RNNdataByPconn.detect_inact_mean = nanmean(-detect_rate);
RNNdataByPconn.detect_inact_std  = nanstd(-detect_rate);

errorbar(cfg.PconnVals,RNNdataByPconn.EA_inact_mean,RNNdataByPconn.EA_inact_std,         ...
         'color',cfg.towersCl,'linewidth',.75)
errorbar(cfg.PconnVals,RNNdataByPconn.detect_inact_mean,RNNdataByPconn.detect_inact_std, ...
         'color',cfg.ctrlCl,'linewidth',.75)
legend({'accum','guided'},'location','best')
set(axs,'xscale','log','xtick',[.05 .5 1])
xlabel('Teacher P_{conn}')
ylabel('Decay rate (%/%)')
xlim([.025 1.5])

%% correlation distributions
fprintf('calculating correlations...\n')

data = load(sprintf('%sdata_and_network_sparse_%g.mat',analysisFilePath,cfg.PconnLevel), ...
              'R2data','Task_data','N','Godata');

% first as is
[EA_corr,detect_corr] = RNNcorrs(data,'RNN');

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')
h(1) = histogram(EA_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(2) = histogram(detect_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(1).EdgeColor = cfg.towersCl;
h(2).EdgeColor = cfg.ctrlCl;
xlabel('Inter-unit r')
ylabel('Prop. unit pairs')
title('RNN')
yl = get(axs, 'ylim');
plot([0 0],[yl(1) yl(2)],'--','color',[.7 .7 .7])
ylim([0 .05])
% legend(h, {'accum.','visually-guided'}, 'location', 'north')

% now add global signal
Go          = data.Godata(1,:,:)+min(data.Godata(:));
Go          = Go./max(Go(:));
data.R2data = data.R2data + cfg.globalSigStrength.*repmat(Go,[size(data.R2data,1) 1 1]);
[EA_corr,detect_corr] = RNNcorrs(data,'RNN');

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')
h(1) = histogram(EA_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(2) = histogram(detect_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(1).EdgeColor = cfg.towersCl;
h(2).EdgeColor = cfg.ctrlCl;
xlabel('Inter-unit r')
ylabel('Prop. unit pairs')
title('RNN + global signal')
yl = get(axs, 'ylim');
plot([0 0],[yl(1) yl(2)],'--','color',[.7 .7 .7])
ylim([0 .3])

% real data
load(summaryFile,'corrSumm')
[EA_corr,detect_corr] = RNNcorrs(corrSumm,'data');

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')
h(1) = histogram(EA_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(2) = histogram(detect_corr(:),cfg.rbins,'DisplayStyle','stairs','Normalization','probability');
h(1).LineWidth = 1;
h(2).LineWidth = 1;
h(1).EdgeColor = cfg.towersCl;
h(2).EdgeColor = cfg.ctrlCl;
xlabel('Inter-ROI r')
ylabel('Prop. ROI pairs')
title('Data')
yl = get(axs, 'ylim');
plot([0 0],[yl(1) yl(2)],'--','color',[.7 .7 .7])
ylim([0 .15])

%% task performance by module (inactivation)
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

load(summaryFile,'RNNdata')
plot([0 inactByConn.percents(end)],[50 50],'--','color',[.7 .7 .7])
h(1) = errorbar(inactByConn.percents,RNNdata.inact_vis.detect_perf,RNNdata.inact_vis.detect_perf_sem,...
               'x','color',cfg.clustCl(1,:),'markerfacecolor',cfg.clustCl(1,:),'markersize',4);
h(2) = errorbar(inactByConn.percents,RNNdata.inact_vis.EA_perf,RNNdata.inact_vis.EA_perf_sem,...
               'o','color',cfg.clustCl(1,:),'markerfacecolor',cfg.clustCl(1,:),'markersize',2.5);
h(3) = errorbar(inactByConn.percents,RNNdata.inact_front.detect_perf,RNNdata.inact_front.detect_perf_sem,...
               'x','color',cfg.clustCl(2,:),'markerfacecolor',cfg.clustCl(2,:),'markersize',4);
h(4) = errorbar(inactByConn.percents,RNNdata.inact_front.EA_perf,RNNdata.inact_front.EA_perf_sem,...
               'o','color',cfg.clustCl(2,:),'markerfacecolor',cfg.clustCl(2,:),'markersize',2.5);
             
thisf = fit(inactByConn.percents,RNNdata.inact_vis.detect_perf',inactByConn.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inactByConn.percents(end),fiteval(thisf,(0:.01:inactByConn.percents(end))'),...
     '--','color',cfg.clustCl(1,:),'linewidth',.75)
thisf = fit(inactByConn.percents,RNNdata.inact_vis.EA_perf',inactByConn.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inactByConn.percents(end),fiteval(thisf,(0:.01:inactByConn.percents(end))'),...
     '-','color',cfg.clustCl(1,:),'linewidth',.75)
thisf = fit(inactByConn.percents,RNNdata.inact_front.detect_perf',inactByConn.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inactByConn.percents(end),fiteval(thisf,(0:.01:inactByConn.percents(end))'),...
     '--','color',cfg.clustCl(2,:),'linewidth',.75)
thisf = fit(inactByConn.percents,RNNdata.inact_front.EA_perf',inactByConn.model, ...
            'startpoint',[50 50 -3],'lower',[40 0 -10],'upper',[100 60 0]);
plot(0:.01:inactByConn.percents(end),fiteval(thisf,(0:.01:inactByConn.percents(end))'),...
     '-','color',cfg.clustCl(2,:),'linewidth',.75)
   
legend(h,{'guided post.','accum. post.','guided front.','accum. front.'});
xlabel('Inactivated units (%)')
ylabel('Perf. (% correct)')
% xlim([0 inact.percents(end)])
xlim([0 10])
ylim([45 100])


%% append bootstrap stats to RNN data
idx = cfg.PconnVals == cfg.PconnLevel;
RNNdata.inact_rand.p_rate_boot = sum(EA_rate(:,idx) > detect_rate(:,idx)) / sum(~isnan(EA_rate(:,idx)));

%% Save data and figure to disk
RNNdataByPconn.xClustCorr_accum_mean = xClustCorr_accum_mean;
RNNdataByPconn.xClustCorr_accum_std  = xClustCorr_accum_mean;
RNNdataByPconn.xClustCorr_guide_mean = xClustCorr_accum_mean;
RNNdataByPconn.xClustCorr_guide_std  = xClustCorr_accum_mean;
RNNdataByPconn.cfg                   = cfg;

save(summaryFile,'RNNdataByPconn','RNNdata','-append')

fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% compute correlations
function [EA_corr,detect_corr] = RNNcorrs(data,datatype)

switch datatype
  case 'RNN'
    pulsesR     = data.R2data(:,:,data.Task_data == 1);
    detectR     = data.R2data(:,:,data.Task_data == 2);
    N           = data.N;

    nPulseTrial = size(pulsesR,3);
    EA_corr     = zeros(N,N,nPulseTrial);
    for iTrial = 1:nPulseTrial
        EA_corr(:,:,iTrial) = corr(pulsesR(:,:,iTrial)');
    end

    nDetectTrial = size(detectR,3);
    detect_corr  = zeros(N,N,nDetectTrial);
    for iTrial = 1:nDetectTrial
        detect_corr(:,:,iTrial) = corr(detectR(:,:,iTrial)');
    end

    EA_corr     = nanmean(EA_corr,3);
    detect_corr = nanmean(detect_corr,3);
    
  case 'data'
    cc    = data.cc_visGuide_wholeTrial_correct;
    nmice = size(cc,3);
    
    for iMouse = 1:nmice
      thiscc = cc(:,:,iMouse); 
      thiscc(logical(eye(size(thiscc,1)))) = nan;
      cc(:,:,iMouse) = thiscc; 
    end
    detect_corr = cc;
    
    cc    = data.cc_accumul_wholeTrial_correct;
    nmice = size(cc,3);
    
    for iMouse = 1:nmice
      thiscc = cc(:,:,iMouse); 
      thiscc(logical(eye(size(thiscc,1)))) = nan;
      cc(:,:,iMouse) = thiscc; 
    end
    EA_corr = cc;
end

end