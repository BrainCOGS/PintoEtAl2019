function taskComp = Pinto2019_figS1_taskPerfComp(analysisFilePath,summaryFile,loadData)

% Pinto2019_figS1_taskPerfComp(analysisFilePath,summaryFile,loadData)
% compares performance between accumulating-towers task and visually guided
% control task
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate all behavioral logs for control trials:
%    concatLogsBehavSummary_owf.m (repo: behaviorAnalysis) (will save owf_concatLog.mat)
%% ------------------------------------------------------------------------

if nargin < 3; loadData = false; end

%% Analysis configuration
cfg.perfThTowers        = 0.6;
cfg.perfThCtrl          = 0.8;
cfg.minNumTrials        = 0;
cfg.minNumTrialsLogReg  = 500;
cfg.towersCl            = widefieldParams.darkgray; 
cfg.ctrlCl              = widefieldParams.darkgreen; 
cfg.memCl               = [149 7 232]./255;
cfg.posBins             = 0:5:300;
cfg.niter               = 1000;
cfg.visGuideExcludeID   = 22; % this mouse was not used in towers
cfg.perfTh              = 0.6;
cfg.nBins               = 4;
cfg.logRegBins          = linspace(10,200,cfg.nBins+1);

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});

if loadData
  load(summaryFile,'logisticRegr','psych','psychByMouse','taskComp');
  
else

  %% get combined log for baseline behavior (opto & wf)
  load([analysisFilePath 'owf_concatLog'],'lg','lg_visGuide','lg_memGuide')

  %% get combined log for baseline behavior (opto & wf)
  lg.firstTrialofBlock          = [true diff(lg.meanPerfBlock)~=0];
  lg_visGuide.firstTrialofBlock = [true diff(lg_visGuide.meanPerfBlock)~=0];
  lg_memGuide.firstTrialofBlock = [true diff(lg_memGuide.meanPerfBlock)~=0];

  %% analyze baseline behavior, towers
  fprintf('calculating psychometrics...\n')
  % psychometric metamouse
  psych        = psychometricFit(lg.choice,lg.nCues_RminusL,true);
  % psychometric avg by mouse
  psychByMouse = xMousePyschometric(lg.choice,lg.nCues_RminusL,lg.mouseID);
  % logistic regression
  fprintf('calculating logistic regression...\n')
  logisticRegr = logRegressionFromConcatLog(lg,cfg.minNumTrialsLogReg,cfg.nBins);

  %% save analysis summary file
  if ~isempty(summaryFile)
    if isempty(dir(summaryFile))
      save(summaryFile,'logisticRegr','psych','psychByMouse','-v7.3')
    else
      save(summaryFile,'logisticRegr','psych','psychByMouse','-append')
    end
  end

  %% continue task comparison
  load(summaryFile, 'taskComp')
  taskComp.motor          = compareTasksMotor(lg,lg_visGuide,lg_memGuide,[],false,false); % motor indicators
  taskComp.motorWF        = compareTasksMotor(lg,lg_visGuide,lg_memGuide,[],false,true); % motor indicators
  taskComp.evidModulation = evModulationTaskComp(lg,lg_visGuide); % modulation of perf by evidence
  save(summaryFile, 'taskComp', '-append')
end

%% Configure figure and panels
layout       = [ 1:4  ...
               ; 5:8  ...
               ; 9:12 ...
               ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% Plot all combined psychometrcis for each mouse
iPanel       = pnlPsychAll(fig, iPanel, psychByMouse, psych);

%% Plot all combined logistic regressions for each mouse
iPanel       = pnlLogRegAll(fig, iPanel, logisticRegr);

%% plot psychometrics  tasks
iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plotPsychometricCurve_ctrl(taskComp.towers.psychPercCorrect, 'bin', axs, cfg.towersCl, true, cfg.towersCl, false);
plotPsychometricCurve_ctrl(taskComp.visGuide.psychPercCorrect, 'bin', axs, cfg.ctrlCl, true, cfg.ctrlCl, false);

[pc,binoerr] = binointervalLP(sum(lg_memGuide.choice == lg_memGuide.trialType),...
                              numel(lg_memGuide.trialType),1-(normcdf(1, 0, 1) - normcdf(-1, 0, 1)));
plot(18,100*pc,'.','color',cfg.memCl,'markersize',7);
errbar(18,100*pc,100*pc-100.*binoerr',cfg.memCl,.75);
set(axs,'xtick',[0:4:16 18],'xticklabel',{'0','','8','','16','memGuide'},'ytick',50:25:100)
ylim([50 100]); xlim([0 19])
ylabel('Perf. (% correct)')

%% plot performance modulation by evidence comparison across mice
iPanel       = pnlTaskComp(fig, iPanel, taskComp, cfg, 'evidModulation'); 

%% plot speed comparison across mice
iPanel       = pnlTaskComp(fig, iPanel, taskComp, cfg, 'speed'); 

%% plot turn onset comparison across mice
iPanel       = pnlTaskComp(fig, iPanel, taskComp, cfg, 'turnOnset'); 

%% plot final decceleration comparison across mice
iPanel       = pnlTaskComp(fig, iPanel, taskComp, cfg, 'accelFinish'); 

%% plot tower-triggered view angle comparison
iPanel       = pnlTowerTrigTheta(fig, iPanel, taskComp, cfg); 

%% plot speed comparison across mice
iPanel       = pnlTaskCompPaired(fig, iPanel, taskComp, cfg, 'speed'); 

%% plot turn onset comparison across mice
iPanel       = pnlTaskCompPaired(fig, iPanel, taskComp, cfg, 'turnOnset'); 

%% plot final decceleration comparison across mice
iPanel       = pnlTaskCompPaired(fig, iPanel, taskComp, cfg, 'accelFinish'); 

%% plot final decceleration comparison across mice
iPanel       = pnlTaskCompPaired(fig, iPanel, taskComp, cfg, 'viewAngSD'); 

%% export
versionInfo.stats = taskComp.stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% -------------------------------------------------------------------
%% psychometrics for all mice on the paper
function iPanel = pnlPsychAll(fig, iPanel, psychByMouse, psych)

% plot best fitting sigmoid for each mouse and overlay metamouse
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);

hold(axs, 'on')

for iMouse = 1:size(psychByMouse.goodCurves,2)
  plot(psychByMouse.fitXaxis,100.*psychByMouse.goodFits(:,iMouse),'-', ...
                'linewidth',         FormatDefaults.linewidthThin,     ...        
                'color'    ,         FormatDefaults.lightGray)
end

plotPsychometricCurve_ctrl(psych, 'all', axs, 'k', true, 'k', false);

set(axs, 'ytick', 0:25:100)
ylabel(axs,'Went right (%)')
xlabel(axs, '\Delta towers (#R - #L)')

end

%% -------------------------------------------------------------------
%% logistic regression for all mice on the paper
function iPanel = pnlLogRegAll(fig, iPanel, logisticRegr)

% plot inidvidual mice and average logistic regression
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);

hold(axs, 'on')

for iMouse = 1:logisticRegr.nGoodMice
    plot(axs, toBinCenters(logisticRegr.bins),                    ...
              logisticRegr.goodValues(:,iMouse),                  ...
              '-','linewidth',  FormatDefaults.linewidthThin,     ...
                  'color',      FormatDefaults.lightGray);
end

errbar(toBinCenters(logisticRegr.bins), logisticRegr.values, logisticRegr.sem, 'k', .75, 0, '-');

xlim([0 200]); ylim([-.1 .4])
set(axs, 'xtick', [0 100 200], 'ytick', 0:.1:.4)
xlabel('Cue y (cm)')
ylabel('Weight on decision')

end

%% -------------------------------------------------------------------
%% task comparison
function iPanel = pnlTaskCompPaired(fig, iPanel, taskComp, cfg, plotWhat)

switch plotWhat
  case 'percCorrect'
    towers   = taskComp.towers.percCorrect.mousePC'; 
    visguide = taskComp.visGuide.percCorrect.mousePC';
    towerAvg = taskComp.towers.percCorrect.mouseAvg;
    towerSem = taskComp.towers.percCorrect.mouseSem;
    vgAvg    = taskComp.visGuide.percCorrect.mouseAvg;
    vgSem    = taskComp.visGuide.percCorrect.mouseSem;
    p        = taskComp.stats.p_percCorrect;
    ylbl     = 'Perf. (% correct)';
    
  case 'evidModulation'
    towers   = taskComp.evidModulation.towers.PC_highMinusLowDelta; 
    visguide = taskComp.evidModulation.visGuide.PC_highMinusLowDelta;
    towerAvg = taskComp.evidModulation.towers.PC_highMinusLowDelta_avg;
    towerSem = taskComp.evidModulation.towers.PC_highMinusLowDelta_sem;
    vgAvg    = taskComp.evidModulation.visGuide.PC_highMinusLowDelta_avg;
    vgSem    = taskComp.evidModulation.visGuide.PC_highMinusLowDelta_sem;
    p        = taskComp.evidModulation.stats.towersVSvisGuide_p;
    ylbl     = sprintf('Perf. (%s %% correct)\n(high - low evidence)','\Delta');
    
  case 'speed'
    towers   = taskComp.motorWF.towers.speed; 
    visguide = taskComp.motorWF.visGuide.speed;
    towerAvg = taskComp.motorWF.towers.speed_avg;
    towerSem = taskComp.motorWF.towers.speed_sem;
    vgAvg    = taskComp.motorWF.visGuide.speed_avg;
    vgSem    = taskComp.motorWF.visGuide.speed_sem;
    p        = taskComp.motorWF.stats.p_speed;
    ylbl     = 'Speed (cm/s)';
    
  case 'turnOnset'
    towers   = taskComp.motorWF.towers.turnOnset; 
    visguide = taskComp.motorWF.visGuide.turnOnset;
    towerAvg = taskComp.motorWF.towers.turnOnset_avg;
    towerSem = taskComp.motorWF.towers.turnOnset_sem;
    vgAvg    = taskComp.motorWF.visGuide.turnOnset_avg;
    vgSem    = taskComp.motorWF.visGuide.turnOnset_sem;
    p        = taskComp.motorWF.stats.p_turnOnset;
    ylbl     = 'Turn onset (cm)';

  case 'accelFinish'
    towers   = taskComp.motorWF.towers.accel_last50; 
    visguide = taskComp.motorWF.visGuide.accel_last50;
    towerAvg = taskComp.motorWF.towers.accel_last50_avg;
    towerSem = taskComp.motorWF.towers.accel_last50_sem;
    vgAvg    = taskComp.motorWF.visGuide.accel_last50_avg;
    vgSem    = taskComp.motorWF.visGuide.accel_last50_sem;
    p        = taskComp.motorWF.stats.p_accel_last50;
    ylbl     = 'Acceleration (last 50 cm)(cm/s^2)';
    
  case 'viewAngSD'
    towers   = taskComp.motorWF.towers.viewAngSD'; 
    visguide = taskComp.motorWF.visGuide.viewAngSD';
    towerAvg = taskComp.motorWF.towers.viewAngSD_avg;
    towerSem = taskComp.motorWF.towers.viewAngSD_sem;
    vgAvg    = taskComp.motorWF.visGuide.viewAngSD_avg;
    vgSem    = taskComp.motorWF.visGuide.viewAngSD_sem;
    p        = taskComp.motorWF.stats.p_viewAngSD;
    ylbl     = 'View Angle S.D. (deg)';
end

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

plot([ones(size(towers)) 1+ones(size(towers))]',[towers visguide]','-','linewidth',.5,'color',[.6 .6 .6])
errbar(0.85,towerAvg,towerSem,cfg.towersCl); 
plot(0.85,towerAvg,'.','color',cfg.towersCl,'markersize',20) 
errbar(2.15,vgAvg,vgSem,cfg.ctrlCl);
plot(2.15,vgAvg,'.','color',cfg.ctrlCl,'markersize',20) 

yl = get(axs,'ylim');
text(1.5,yl(2)*.95,sprintf('p = %1.2g',p),'color','k','horizontalAlignment','center','fontsize',9)

set(axs,'xtick',1:2,'xticklabel',{'accum.-towers';'vis.-guided'})
rotateXLabels(axs,45)
ylabel(ylbl)
xlim([.65 2.35])

end

%% -------------------------------------------------------------------
%% task comparison
function iPanel = pnlTaskComp(fig, iPanel, taskComp, cfg, plotWhat)

switch plotWhat
  case 'percCorrect'
    towers   = taskComp.towers.percCorrect.mousePC'; 
    visguide = taskComp.visGuide.percCorrect.mousePC';
    memguide = taskComp.memGuide.percCorrect.mousePC';
    p        = taskComp.stats.p_percCorrect;
    ylbl     = 'Perf. (% correct)';
    
  case 'evidModulation'
    towers   = taskComp.evidModulation.towers.PC_highMinusLowDelta; 
    visguide = taskComp.evidModulation.visGuide.PC_highMinusLowDelta;
    p        = taskComp.evidModulation.stats.towersVSvisGuide_p;
    ylbl     = sprintf('Perf. (%s %% correct)\n(high - low evidence)','\Delta');
    
  case 'speed'
    towers   = taskComp.motor.towers.speed; 
    visguide = taskComp.motor.visGuide.speed;
    memguide = taskComp.motor.memGuide.speed;
    p        = taskComp.motor.stats.p_speed;
    ylbl     = 'Speed (cm/s)';
    
  case 'turnOnset'
    towers   = taskComp.motor.towers.turnOnset; 
    visguide = taskComp.motor.visGuide.turnOnset;
    memguide = taskComp.motor.memGuide.turnOnset;
    p        = taskComp.motor.stats.p_turnOnset;
    ylbl     = 'Turn onset (cm)';

  case 'accelFinish'
    towers   = taskComp.motor.towers.accel_last50; 
    visguide = taskComp.motor.visGuide.accel_last50;
    memguide = taskComp.motor.memGuide.accel_last50;
    p        = taskComp.motor.stats.p_accel_last50;
    ylbl     = 'Acceleration (last 50 cm)(cm/s^2)';
end

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

switch plotWhat
  case 'evidModulation'
    boxplot([visguide; towers],[ones(size(visguide)); 1+ones(size(towers))],...
            'colors',[cfg.ctrlCl; cfg.towersCl])
    yl = get(axs,'ylim');
    text(1.5,yl(2)*.95,sprintf('p = %1.2g',p),'color','k','horizontalAlignment','center','fontsize',9)
    set(axs,'xtick',1:3,'xticklabel',{'vis.-guided','accum.-towers'})
    rotateXLabels(axs,45)
    ylabel(ylbl)
    
  otherwise
    boxplot([visguide; memguide ; towers],[ones(size(visguide)); 1+ones(size(memguide)); 2+ones(size(towers))],...
            'colors',[cfg.ctrlCl; cfg.memCl; cfg.towersCl])
    yl = get(axs,'ylim');
    text(2,yl(2)*.95,sprintf('p = %1.2g',p),'color','k','horizontalAlignment','center','fontsize',9)
    set(axs,'xtick',1:3,'xticklabel',{'vis.-guided','mem.-guided','accum.-towers'})
    rotateXLabels(axs,45)
    ylabel(ylbl)
end
box off
end

%% -------------------------------------------------------------------
%% compare tower-triggered view angle
function iPanel       = pnlTowerTrigTheta(fig, iPanel, taskComp, cfg)

va_towers_avg  = taskComp.motorWF.towers.towerTrigVAng.alignedVAng_sideAvg_avg;
va_towers_sem  = taskComp.motorWF.towers.towerTrigVAng.alignedVAng_sideAvg_sem;
va_vg_avg      = taskComp.motorWF.visGuide.towerTrigVAng.alignedVAng_sideAvg_avg;
va_vg_sem      = taskComp.motorWF.visGuide.towerTrigVAng.alignedVAng_sideAvg_sem;
xaxis          = taskComp.motorWF.towers.towerTrigVAng.posAxis;

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
hold(axs, 'on')

%%
% figure; hold on; axs = gca;
plot([xaxis(1) xaxis(end)],[0 0],'--','linewidth',.5,'color',[.7 .7 .7])

plot(xaxis,va_towers_avg,'-','linewidth',1,'color',cfg.towersCl) 
plot(xaxis,va_towers_avg-va_towers_sem,'--','linewidth',.25,'color',cfg.towersCl) 
plot(xaxis,va_towers_avg+va_towers_sem,'--','linewidth',.25,'color',cfg.towersCl) 

plot(xaxis,va_vg_avg,'-','linewidth',1,'color',cfg.ctrlCl) 
plot(xaxis,va_vg_avg-va_vg_sem,'--','linewidth',.25,'color',cfg.ctrlCl) 
plot(xaxis,va_vg_avg+va_vg_sem,'--','linewidth',.25,'color',cfg.ctrlCl) 

ylabel('\Theta (deg)')
xlabel('y pos from tower (cm)')

yl = get(axs,'ylim');
ptask = taskComp.motorWF.stats.towerTrigVAng_ANOVA_p(1);
ppos  = taskComp.motorWF.stats.towerTrigVAng_ANOVA_p(2);
text(xaxis(1) + 5,yl(2)*.95,sprintf('%s = %1.2g\n%s = %1.2g','p_{task}',ptask,'p_{pos}',ppos), ...
  'color','k','horizontalAlignment','left','fontsize',9)
xlim([xaxis(1) xaxis(end)])

end

%% -------------------------------------------------------------------
%% compare difference in performance between high and low evidence trials
function evidModulation = evModulationTaskComp(lg,lg_visGuide)

% vis guide
mids       = unique(lg_visGuide.mouseID);
lowDelta   = zeros(numel(mids),1);
highDelta  = lowDelta;

for iMouse = 1:numel(mids)
  [ch,tt,delta]     = selectMouseTrials(lg_visGuide, mids(iMouse), 'choice','trialType','nCues_RminusL');
  delta             = abs(delta);
  lowDelta(iMouse)  = 100*(sum(ch(delta<median(delta)) == tt(delta<median(delta)))/sum(delta<median(delta)));
  highDelta(iMouse) = 100*(sum(ch(delta>median(delta)) == tt(delta>median(delta)))/sum(delta>median(delta)));
end

evidModulation.visGuide.PC_highMinusLowDelta     = highDelta - lowDelta;
evidModulation.visGuide.PC_highMinusLowDelta_avg = mean(highDelta - lowDelta);
evidModulation.visGuide.PC_highMinusLowDelta_sem = std(highDelta - lowDelta) ./ sqrt(numel(mids)-1);

if lillietest(highDelta - lowDelta)
  evidModulation.visGuide.PC_highMinusLowDelta_test  = 'signrank';
  evidModulation.visGuide.PC_highMinusLowDelta_p     = signrank(highDelta - lowDelta);
else
  evidModulation.visGuide.PC_highMinusLowDelta_test  = 'ttest';
  [~,evidModulation.visGuide.PC_highMinusLowDelta_p] = ttest(highDelta - lowDelta);
end

% towers
mids       = unique(lg.mouseID);
lowDelta   = zeros(numel(mids),1);
highDelta  = lowDelta;

for iMouse = 1:numel(mids)
  [ch,tt,delta]     = selectMouseTrials(lg, mids(iMouse), 'choice','trialType','nCues_RminusL');
  delta             = abs(delta);
  lowDelta(iMouse)  = 100*(sum(ch(delta<median(delta)) == tt(delta<median(delta)))/sum(delta<median(delta)));
  highDelta(iMouse) = 100*(sum(ch(delta>median(delta)) == tt(delta>median(delta)))/sum(delta>median(delta)));
end

evidModulation.towers.PC_highMinusLowDelta     = highDelta - lowDelta;
evidModulation.towers.PC_highMinusLowDelta_avg = mean(highDelta - lowDelta);
evidModulation.towers.PC_highMinusLowDelta_sem = std(highDelta - lowDelta) ./ sqrt(numel(mids)-1);

if lillietest(highDelta - lowDelta)
  evidModulation.towers.PC_highMinusLowDelta_test  = 'signrank';
  evidModulation.towers.PC_highMinusLowDelta_p     = signrank(highDelta - lowDelta);
else
  evidModulation.towers.PC_highMinusLowDelta_test  = 'ttest';
  [~,evidModulation.towers.PC_highMinusLowDelta_p] = ttest(highDelta - lowDelta);
end

% task comp
if lillietest(evidModulation.towers.PC_highMinusLowDelta) || lillietest(evidModulation.visGuide.PC_highMinusLowDelta)
  evidModulation.stats.towersVSvisGuide_test  = 'signrank';
  evidModulation.stats.towersVSvisGuide_p     = ranksum(evidModulation.visGuide.PC_highMinusLowDelta,evidModulation.towers.PC_highMinusLowDelta);
else
  evidModulation.stats.towersVSvisGuide_test  = 'ttest2';
  [~,evidModulation.stats.towersVSvisGuide_p] = ttest2(evidModulation.visGuide.PC_highMinusLowDelta,evidModulation.towers.PC_highMinusLowDelta);
end

end