function Pinto2019_figS1_speedSingleTrialEg(analysisFilePath)

% Pinto2019_figS1_speedSingleTrialEg(analysisFilePath)
% plots single-trial examples of speed and view angle for tasks
% analysisFilePath is path for data analysis files to be loaded

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To generate all behavioral logs for control trials:
%    concatLogsBehavSummary_owf.m (repo: behaviorAnalysis) (will save owf_concatLog.mat)
%% ------------------------------------------------------------------------
  
%% Analysis configuration
cfg.towersCl            = widefieldParams.darkgray; %[0 0 0];
cfg.ctrlCl              = widefieldParams.darkgreen; %[.6 .6 .6];
cfg.memCl               = [149 7 232]./255;

cfg.lStart             = 30;
cfg.lCue               = 200;
cfg.lMemory            = 100;
cfg.lMaze              = cfg.lCue + cfg.lMemory;
cfg.posBins            = 0:5:cfg.lMaze;

cfg.egMouseID          = 45;
cfg.numTrialsEg        = 20;
cfg.doSmooth           = 1; % smooth speed traces for display?
cfg.smoothWin          = 5; % time points
cfg.firstTrial         = [105 7 1]; %
cfg.choiceSide         = 0;
cfg.taskList           = {'towers','ctrl','mem'};

%% Version tracking
fprintf('collecting version info...\n')
codeFile         = mfilename('fullpath');
versionInfo      = collectVersionInfo(codeFile, cfg, [], [], {});


%% get combined log for baseline behavior (opto & wf)
load([analysisFilePath 'owf_concatLog'],'lg','lg_visGuide')

%% Configure figure and panels
layout       = 1:4 ;
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

% Plot single trial traces
iPanel       = pnlEgSession(fig, iPanel, lg, cfg, 'towers', 'speed');
iPanel       = pnlEgSession(fig, iPanel, lg_visGuide, cfg, 'ctrl', 'speed');
iPanel       = pnlEgSession(fig, iPanel, lg, cfg, 'towers', 'viewAngle');
iPanel       = pnlEgSession(fig, iPanel, lg_visGuide, cfg, 'ctrl', 'viewAngle');

%% export
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])
end

%% --------------------------------------------------------------------
%% Plot speed or view angle traces
function iPanel = pnlEgSession(fig, iPanel, lg, cfg, whichTask, whichTrace)

%%
switch whichTask
  case 'towers'
    lbl     = 'Accum.-towers task';
  case 'ctrl'
    lbl     = 'Vis.-guided task';
  case 'mem'
    lbl     = 'Mem.-guided task';
end

idx1    = strcmp(cfg.taskList,whichTask);

pos     = lg.pos(lg.mouseID == cfg.egMouseID & lg.choice == cfg.choiceSide & lg.choice == lg.trialType & lg.excessTravel < .1);
pos     = pos(cfg.firstTrial(idx1):cfg.numTrialsEg+cfg.firstTrial(idx1)-1);

switch whichTrace
  case 'speed' 
    time    = lg.time(lg.mouseID == cfg.egMouseID & lg.choice == cfg.choiceSide & lg.choice == lg.trialType & lg.excessTravel < .1);
    time    = time(cfg.firstTrial(idx1):cfg.numTrialsEg+cfg.firstTrial(idx1)-1);
    for iT = 1:numel(time)
      displ        = [0 0; diff(pos{iT}(:,1:2))];
      pos{iT}(:,3) = sqrt(sum(displ.^2,2)) ./ [0; diff(time{iT}(1:size(displ,1)))]; 
    end
    ylbl    = 'Speed (cm/s)';
    yl      = [0 125];
  case 'viewAngle'
    ylbl    = 'View angle (deg)';
    if cfg.choiceSide == 0
      yl      = [-100 25];
    else
      yl      = [-25 100];
    end
end

data        = sampleViewAngleVsY(pos, cfg.posBins, false);

%smooth?
if cfg.doSmooth
  dat = cell(1,size(data,2));
  for iT = 1:size(data,2); dat{iT} = data(:,iT); end
  w    = num2cell(ones(1,numel(dat)).*cfg.smoothWin);
  data = cell2mat((cellfun(@smooth,dat,w,'UniformOutput',false)));
end

iPanel      = iPanel + 1;
axs         = fig.panel(iPanel);
cl = cfg.([whichTask 'Cl']);
plot(cfg.posBins,data','color',cl,'linewidth',.25)


% axis   tight
ylim(yl)
xlim([0 300])
xlabel (axs, 'y pos (cm)')
ylabel (axs, ylbl)
title(lbl)
            
end