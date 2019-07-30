function Pinto2019_fig3_WFdynamics(analysisFilePath,summaryFile)

% Pinto2019_fig3_WFdynamics(analysisFilePath,summaryFile)
% plots avg activity dynamics results for widefield data, both tasks
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 0) A list of all the used recordings (relative paths) can be found in the
%    class object widfield_recLs.m (repo: imaging)
% 1) To preprocess widefield data (in this sequence): 
%    1.1) widefieldPreProc.m (to generate dff.mat, info.mat, repo: imaging)
%    1.2) summarizeVirmenLog_widefield.m (to sync with behav - logSumm.mat, repo: widfieldImaging)
%    1.3) to get ROIs, the following must be done once for each mouse: (repo: widfieldImaging)
%       - choose a reference image and save it as variable name meanproj in
%         a file named refIm.mat, saved to that mouse's root directory
%       - use setBregmaAndMidline to set bregma and midline location,
%         append output to refIm.mat
%       - analyze maps: airpuffsMap.m & visualCortexMap.m
%       - register maps to each other: summarizeSensoryMaps.m
%       - register sensory maps and Allen Brain Atlas: registerMapsToAllen.m
%         to obtain a refernce ROI map, saved to refROI.mat
%    1.4) to generate ROI dff, extractDffFromROI.m
% 2) To generate trial-triggered averages and general stats: (repo: widfieldImaging)
%    2.1) avgDffByTrialEpoch.m for pixel trial averages
%    2.2) avgDffROI.m for ROI trial averages
%    2.3) dffANOVA3widefield.m for task modulation analysis
% 3) To summarize experiments:
%    summarizeAvgDff.m (will save avgDffSummary.mat, loaded here)(repo: widfieldImaging)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.egRec            = 'ai3/20170126';
cfg.egRecRawTrace    = 'ai2/20170126';
cfg.egRecFrame1      = 101; % for raw traces with no behav
cfg.egRecDurMin      = 3;  % for raw traces with no behav
cfg.egAddBehav       = true; % to also plot behav events
cfg.egNtrials        = 10;
cfg.egFirstTrialID   = 98;
cfg.egROI            = {'VISp-R','VISa-R','RSP-R','mMOs-R','MOp-R'};
cfg.egTrialType      = 'correctL';
cfg.ROIcontraOnly    = true; % to plot just contralateral ROIs in heatmap
cfg.frameBinsEg      = 0:50:300;
cfg.posBins          = 0:5:300;
cfg.dffCMap          = parula;
cfg.normROImaps      = true; % normalize activity in ROI avg plots
cfg.compressPctile   = []; % compress dff scale in pxl maps for display
cfg.normPxlMaps      = true; % normalize activity?
cfg.ROI_avgOReg      = 'lines'; % what type of plot
cfg.alpha            = 0.01; % FDR significance threshold

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/avgDffSummary.mat'],'avgDffSumm');

%% compile stats 
stats.avgDff  = avgDffSumm.stats;
stats.avgDff  = taskANOVA(avgDffSumm,cfg);

%% retrieve egs
egdata        = load([analysisFilePath cfg.egRecRawTrace '/dffROI.mat']);
egavg         = load([analysisFilePath cfg.egRec '/avgDff.mat']);
if strcmpi(cfg.ROI_avgOReg,'eg')
  egavgROI    = load([analysisFilePath cfg.egRec '/avgROI.mat']);
else
  egavgROI    = [];
end
if cfg.egAddBehav; egLg = load([analysisFilePath cfg.egRecRawTrace '/behavLog.mat']); end

%% Configure figure and panels
layout       = [ 1  3  3  3  3  3  4  ...
               ; 2  3  3  3  3  3  4  ...
               ; 5  6  7  8  9  10 11 ...
               ; 12 13 14 15 16 17 18 ...
               ; 19 20 20 21 22 22 23 ...
               ; 24 25 26 27 28 28 28 ...
               ];
fig          = PaneledFigure(layout, 'tiny');
iPanel       = 0;

%% Plot rig schematics
iPanel          = pnlRigScheme(fig, iPanel, cfg);
iPanel          = iPanel+1; 

%% Plot eg ROI traces
fig3stats = [];
if cfg.egAddBehav
  [iPanel,fig3stats]        = pnlRawTracesWithBehav(fig, iPanel, egdata, egLg.logSumm, cfg, fig3stats);
else
  iPanel        = pnlRawTraces(fig, iPanel, egdata, cfg);
  iPanel        = iPanel + 1;
end
%% plot towers task maze schematic
iPanel          = pnlMzScheme(fig, iPanel, cfg, 'towers');

%% time progression, towers task
% passing the data from the other maze will result in joint normalization
[iPanel,fig3stats]   = pnlFramesEg(fig, iPanel, egavg.avgDff.maze(end), cfg, [], fig3stats,'towers');

%% plot towers task maze schematic
iPanel          = pnlMzScheme(fig, iPanel, cfg, 'visualGuide');

%% time progression, visually guided task
[iPanel,fig3stats]  = pnlFramesEg(fig, iPanel, egavg.avgDff.maze(1), cfg, [], fig3stats, 'visGuide');

%% plot towers task maze schematic
iPanel          = pnlMzScheme(fig, iPanel, cfg, 'towers');

%% time progression, towers task, ROI
switch cfg.ROI_avgOReg
  case 'eg'
    iPanel      = pnlFramesROIEg(fig, iPanel, egavgROI.avgROI.maze(end), cfg, avgDffSumm.ROIlbl, sortidx);
  case 'avg'
    iPanel      = pnlFramesROIAvg(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, sortidx, 'accumul');
  case 'lines'
    [iPanel,fig3stats]  = pnlFramesROIAvgLines(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, [], 'accumul', stats.avgDff, fig3stats );
end

%% plot towers task maze schematic
iPanel          = pnlMzScheme(fig, iPanel, cfg, 'visualGuide');

%% time progression, control task, ROI
switch cfg.ROI_avgOReg
  case 'eg'
    iPanel      = pnlFramesROIEg(fig, iPanel, egavgROI.avgROI.maze(1), cfg, avgDffSumm.ROIlbl, sortidx);
  case 'avg'
    iPanel      = pnlFramesROIAvg(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, sortidx, 'visGuide');
  case 'lines'
    [iPanel,fig3stats]   = pnlFramesROIAvgLines(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, [], 'visGuide', stats.avgDff, fig3stats );
end% 

%% insets
iPanel          = iPanel+2;
iPanel          = pnlFramesROIAvgLinesInset(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, [], 'accumul');
iPanel          = pnlFramesROIAvgLinesInset(fig, iPanel, avgDffSumm, cfg, avgDffSumm.ROIlbl, [], 'visGuide');

%% (strip down a bit for saving)
WFstats         = stats;
avgDffSumm      = rmfield(avgDffSumm,'mouse');

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'avgDffSumm','WFstats','fig3stats','-v7.3')
  else
    save(summaryFile,'avgDffSumm','WFstats','fig3stats','-append')
  end
end

%% Save figure to disk
versionInfo.stats = stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% rig schematics
function iPanel = pnlRigScheme(fig,iPanel,cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/rigWf.png']);
imshow(im); axis off; axis image

end

%% ------------------------------------------------------------------------
%% raw traces
function iPanel  = pnlRawTraces(fig, iPanel, data, cfg)

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

ROI       = cellfun(@(x)(find(strcmpi(data.ROIlbl,x))),cfg.egROI);
frameRate = 10;
frames    = cfg.egRecFrame1:cfg.egRecFrame1+cfg.egRecDurMin*60*frameRate;
offset    = .07;

for iROI  = 1:numel(ROI)
  [cl,lbl] = getDefaultROIcl(cfg.egROI{iROI});
  plot(data.dffROI(frames,ROI(iROI))+(iROI-1)*offset,'-','color',cl,'linewidth',1)
  text(numel(frames)*.88,(iROI-1)*offset+.03,lbl,'color',cl,'fontsize',9)
end

set(axs,'xtick',[],'xcolor','w','ytick',[],'ycolor','w')

% time scale, dff scale
axis tight
yl = get(axs,'ylim');
plot([1 1+20*frameRate],[yl(1) yl(1)],'k-')
text(mean([1 1+20*frameRate]), yl(1)-.01, ...
     '20 s', 'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 10)
plot([1 1],[yl(1) yl(1)+.05],'k-')
th = text(-20, yl(1)+.025, '\DeltaF/F (0.05)', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 10);
set(th,'rotation',90)
        
end

%% ------------------------------------------------------------------------
%% raw traces with behavioral events
function [iPanel,figstats]  = pnlRawTracesWithBehav(fig, iPanel, data, lg, cfg, figstats)

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

% get ROI dff
ROI        = cellfun(@(x)(find(strcmpi(data.ROIlbl,x))),cfg.egROI);
frameRate  = 10;
offset     = .07;
frames     = (lg.binned.camFrameID{cfg.egFirstTrialID}(1):lg.binned.camFrameID{cfg.egFirstTrialID+cfg.egNtrials-1}(end))';
scale      = (max(data.dffROI(frames,ROI(1))) - min(data.dffROI(frames,ROI(1))))/1.5;

% make behavioral event vectors
ypos       = zeros(size(frames));
viewang    = zeros(size(frames));
speed      = zeros(size(frames));
trialStart = zeros(size(frames));
cueStart   = zeros(size(frames));
delayStart = zeros(size(frames));
rw         = zeros(size(frames));
timeout    = zeros(size(frames));
towersR    = zeros(size(frames));
towersL    = zeros(size(frames));

for iTrial = cfg.egFirstTrialID:cfg.egFirstTrialID+cfg.egNtrials-1
  frameID             = lg.binned.camFrameID{iTrial} - frames(1) + 1;
  thisy               = lg.binned.pos{iTrial}(:,2);
  ypos(frameID(1:numel(thisy)))       = thisy;
  viewang(frameID(1:numel(thisy)))    = lg.binned.pos{iTrial}(:,3);
  thisdispl                           = sqrt(lg.binned.pos{iTrial}(:,1).^2 + lg.binned.pos{iTrial}(:,2).^2);
  thisspeed                           = abs([0; diff(thisdispl)] ./ [0; diff(lg.binned.time{iTrial}(1:numel(thisy)))]);
  thisspeed(isnan(thisspeed))         = 0;
  speed(frameID(1:numel(thisy)))      = smooth(thisspeed);
  startFID            = lg.camFrameNum{iTrial}(find(lg.pos{iTrial}(:,2)>=-30,1,'first')) - frames(1) + 1;
  cueFID              = lg.camFrameNum{iTrial}(lg.keyFrames{iTrial}(1)) - frames(1) + 1;
  delFID              = lg.camFrameNum{iTrial}(lg.keyFrames{iTrial}(2)) - frames(1) + 1;
  rwFID               = lg.camFrameNum{iTrial}(lg.keyFrames{iTrial}(4)) - frames(1) + 1;
  cueStart(cueFID)    = 1;
  delayStart(delFID)  = 1;
  trialStart(startFID)= 1;
  if lg.trialType(iTrial) == lg.choice(iTrial)
    rw(rwFID)         = 1;
  else
    timeout(rwFID)    = 1;
  end
  towerFID            = zeros(1,numel(lg.cueOnset_R{iTrial}));
  for iT = 1:numel(lg.cueOnset_R{iTrial})
    thist             = lg.cueOnset_R{iTrial}(iT);
    towerFID(iT)      = lg.camFrameNum{iTrial}(find(lg.time{iTrial}>=thist,1,'first')) - frames(1) + 1;
  end
  towersR(towerFID)   = 1;
  towerFID            = zeros(1,numel(lg.cueOnset_L{iTrial}));
  for iT = 1:numel(lg.cueOnset_L{iTrial})
    thist             = lg.cueOnset_L{iTrial}(iT);
    towerFID(iT)      = lg.camFrameNum{iTrial}(find(lg.time{iTrial}>=thist,1,'first')) - frames(1) + 1;
  end
  towersL(towerFID)   = 1;
end

% plot behavioral events at bottom: speed, view angle, y pos
figstats.panelB.ypos        = ypos;
figstats.panelB.viewAngle   = viewang;
figstats.panelB.speed       = speed;
figstats.panelB.rightTowers = towersR;
figstats.panelB.leftTowers  = towersL;
figstats.panelB.reward      = rw;
figstats.panelB.timeout     = timeout;
figstats.panelB.xaxis       = 0:1/frameRate:numel(frames)/frameRate-1/frameRate;

maxy     = max(ypos);
ypos     = ypos ./ maxy .* scale/2;
maxva    = max(abs(viewang));
viewang  = viewang ./ maxva .* scale/2;
maxspeed = max(speed);
speed    = speed ./ maxspeed .* scale/2;

plot(ypos,'-','color',[1 1 1]./2)
text(numel(frames)+30,0,'y pos','color',[1 1 1]./2,'fontsize',9);
plot(speed+offset*.75,'-','color',[1 1 1]./2)
text(numel(frames)+30,offset*.8,'speed','color',[1 1 1]./2,'fontsize',9);
plot(viewang+offset*1.5,'-','color',[1 1 1]./2)
plot(zeros(size(viewang))+offset*1.5,':','color',[1 1 1] - .3)
text(numel(frames)+30,offset*1.6,'view angle (\theta)','color',[1 1 1]./2,'fontsize',9);

% plot ROIs
for iROI  = 1:numel(ROI)
  [cl,lbl] = getDefaultROIcl(cfg.egROI{iROI});
  plot(data.dffROI(frames,ROI(iROI))+(iROI-1)*offset+offset*2.5,'-','color',cl,'linewidth',1)
  text(numel(frames)*.03,(iROI-1)*offset-.03+offset*2.5,lbl,'color',cl,'fontsize',9)
  figstats.panelB.(['dff_' cfg.egROI{iROI}(1:end-2)]) = data.dffROI(frames,ROI(iROI));
end

% plot reward, timeout, cue, delay and towers
yl = get(axs,'ylim');
tR = find(towersR==1);
tL = find(towersL==1);
for iR = 1:numel(tR)
  plot(tR(iR),yl(2)*.97,'>','color',widefieldParams.myblue,'markersize',2.5,'linewidth',.25) %'markerfacecolor',widefieldParams.myblue,
end
for iL = 1:numel(tL)
  plot(tL(iL),yl(2),'<','color',widefieldParams.myred,'markersize',2.5,'linewidth',.25) % 'markerfacecolor',widefieldParams.myred,
end

rwt = find(rw==1);
for iR = 1:numel(rwt)
  text(rwt(iR),yl(2)+.01,'\downarrow','color','k','fontsize',12,'fontweight','bold')
end
tot = find(timeout==1);
for iR = 1:numel(tot)
  text(tot(iR),yl(2)+.01,'\downarrow','color',FormatDefaults.errorCl,'fontsize',12,'fontweight','bold')
end

set(axs,'xtick',[],'xcolor','w','ytick',[],'ycolor','w')

% time scale, dff scale, behavioral variable scales
axis tight
yl = get(axs,'ylim');
plot([-15 20*frameRate-15],[yl(1)-.01 yl(1)-.01],'k-')
text(mean([1 1+20*frameRate]), yl(1)-.04, ...
     '20 s', 'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8)
plot([-15 -15],[yl(1)+offset*2.5 yl(1)+.05+offset*2.5],'k-')
th = text(-45, yl(1)+.025+offset*2.5, '\DeltaF/F (0.05)', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8);
set(th,'rotation',90)
plot([-15 -15],[0; max(ypos)/3],'k-')
th = text(-45, max(ypos)/3/2, '1 m', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8);
set(th,'rotation',90)
plot([-15 -15],[offset*.75; offset*.75+20*scale/maxspeed],'k-')
th = text(-45, offset*.75+20*scale/2/maxspeed/2, '20 cm/s', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8);
set(th,'rotation',90)
plot([-15 -15],[offset*1.5; offset*1.5+50*scale/maxva],'k-')
th = text(-45, offset*1.5+20*scale/maxva/2, '50 deg', ...
          'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 8);
set(th,'rotation',90)

%% plot caption
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')
axis off
xlim([0 1]); ylim([0 1]); 
plot(.1,.9,'>','color',widefieldParams.myblue,'markersize',5)
text(.28,.9,'right tower','color',widefieldParams.myblue,'fontsize',10)
plot(.1,.8,'<','color',widefieldParams.myred,'markersize',5)
text(.28,.8,'left tower','color',widefieldParams.myred,'fontsize',10)
text(.1,.7,'\downarrow     reward','color','k','fontsize',10)
text(.1,.6,'\downarrow     time out','color',FormatDefaults.errorCl,'fontsize',10)
% text(.1,.5,'\downarrow     trial start','color',[.6 .6 .6],'fontsize',10)

end

%% ------------------------------------------------------------------------
%% maze schematics
function iPanel = pnlMzScheme(fig,iPanel,cfg,whichOne)

switch whichOne
  case 'towers'
    im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
  case 'visualGuide'
    im     = imread([cfg.stockFigPath '/mazeSchematicVisualGuide.png']);
end

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
imshow(im); axis off; axis image

end

%% ------------------------------------------------------------------------
%% pxl frames
function [iPanel,figstats] = pnlFramesEg(fig,iPanel,data,cfg,data2,figstats,whichData)

if nargin < 5; data2=[]; end

trialType  = cfg.egTrialType;
tBins      = cfg.frameBinsEg;
cmap       = cfg.dffCMap;
cmap(1,:)  = [0 0 0];
  
dff        = data.(trialType).space.frames;
[nX,nY,nZ] = size(dff);

if ~isempty(data2)
  dff2     = data2.(trialType).space.frames;
else
  dff2     = [];
end

if cfg.normPxlMaps
  dff  = cat(3,dff,dff2);
  nZ2  = size(dff,3);
  dff  = dff - repmat(min(dff,[],3),[1 1 nZ2]);
  dff  = dff ./ repmat(max(dff,[],3),[1 1 nZ2]);
  if ~isempty(dff2); dff2 = dff(:,:,nZ+1:nZ2); end
  dff  = dff(:,:,1:nZ);
end

newdff     = zeros(nX,nY,numel(tBins)-1);
taxis      = linspace(0,300,nZ);
for iBin = 1:numel(tBins)-1
  idx              = taxis >= tBins(iBin) & taxis < tBins(iBin+1);
  newdff(:,:,iBin) = nanmean(dff(:,:,idx),3);
end

if ~isempty(dff2)
  newdff2    = zeros(nX,nY,numel(tBins)-1);
  for iBin = 1:numel(tBins)-1
    idx               = taxis >= tBins(iBin) & taxis < tBins(iBin+1);
    newdff2(:,:,iBin) = nanmean(dff2(:,:,idx),3);
  end
end

if ~isempty(cfg.compressPctile)
  cap                     = prctile(newdff(:),cfg.compressPctile);
  newdff(newdff < cap(1)) = cap(1);
  newdff(newdff > cap(2)) = cap(2);
end

iPanel          = iPanel + 1;
for iFrame = 1:size(newdff,3)
  axs    = fig.panel(iPanel);
  if cfg.normPxlMaps
    imagesc(newdff(:,:,iFrame),[-.001 1]);
  else
    imagesc(newdff(:,:,iFrame),[min(newdff(:))-.0001 max(newdff(:))]);
  end
  colormap(axs,cmap); axis image; axis off
  text(5,-20,sprintf('%d-%d cm',tBins(iFrame),tBins(iFrame+1)),'color','k','fontsize',10)
  
  iPanel = iPanel + 1;
end

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\DeltaF/F';

iPanel            = iPanel - 1;

switch whichData
  case 'visGuide'
    figstats.panelC.frames    = newdff;
    figstats.panelC.posBins   = tBins;
    figstats.panelC.taskLabel = 'Vis.-guided';
  case 'towers'
    figstats.panelD.frames  = newdff;
    figstats.panelD.posBins = tBins;
    figstats.panelD.taskLabel = 'Accum.-towers';
end
end

%% ------------------------------------------------------------------------
%% ROI frames
function iPanel = pnlFramesROIEg(fig,iPanel,data,cfg,ROIlbl,sortidx)

if nargin < 6; sortidx = []; end

trialType  = cfg.egTrialType;
taxis      = cfg.posBins;
cmap       = cfg.dffCMap;  
dff        = data.(trialType).space.frames;
[~,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.ROIcontraOnly 
  side   = cfg.egTrialType(end);
  idx    = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl = ROIlbl(idx);
  dff    = dff(:,idx);
end

if cfg.normROImaps
  nt       = size(dff,1);
  minv     = min(dff);
  dff      = (dff - repmat(minv,[nt 1]));
  maxv     = max(dff);
  dff      = dff ./ repmat(maxv,[nt 1]);
end

if ~isempty(sortidx)
  dff    = dff(:,sortidx);
  ROIlbl = ROIlbl(sortidx);
end

iPanel     = iPanel + 1;
axs        = fig.panel(iPanel);
imagesc(taxis,1:numel(ROIlbl),dff');
colormap(axs,cmap); 
set(axs,'ytick',1:numel(ROIlbl),'yticklabel',ROIlbl,'xtick',0:50:300)
axis xy
xlabel('y pos (com)')
 
cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = 'Norm. \DeltaF/F';
% ch.Tick         = [min(dff(:)) max(dff(:))];%[prctile(newdff(:),1) prctile(newdff(:),99)];

end

%% avg ROI activity
function iPanel = pnlFramesROIAvg(fig,iPanel,data,cfg,ROIlbl,sortidx,whichMaze)

if nargin < 6; sortidx   = [];        end
if nargin < 7; whichMaze = 'accumul'; end

trialType  = cfg.egTrialType;
taxis      = cfg.posBins;
cmap       = cfg.dffCMap;  
dff        = data.(['ROIavg_' whichMaze '_' trialType]); 
[~,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.ROIcontraOnly 
  side   = cfg.egTrialType(end);
  idx    = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl = ROIlbl(idx);
  dff    = dff(:,idx,:);
end

if cfg.normROImaps
  nt     = size(dff,1);
  minv   = min(dff);
  dff    = (dff - repmat(minv,[nt 1 1]));
  maxv   = max(dff);
  dff    = dff ./ repmat(maxv,[nt 1 1]);
end

dff = mean(dff,3);

if ~isempty(sortidx)
  dff    = dff(:,sortidx);
  ROIlbl = ROIlbl(sortidx);
end

iPanel     = iPanel + 1;
axs        = fig.panel(iPanel);
if cfg.normPxlMaps
  imagesc(taxis,1:numel(ROIlbl),dff',[0 1]);
else
  imagesc(taxis,1:numel(ROIlbl),dff');
end
colormap(axs,cmap); 
set(axs,'ytick',1:numel(ROIlbl),'yticklabel',ROIlbl,'xtick',0:50:300)
axis xy
xlabel('y pos (cm)')
 
cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = 'Norm. \DeltaF/F';
% ch.Tick         = [min(dff(:)) max(dff(:))];%[prctile(newdff(:),1) prctile(newdff(:),99)];

end

%% avg ROI activity, lines
function [iPanel,figstats] = pnlFramesROIAvgLines(fig,iPanel,data,cfg,ROIlbl,sortidx,whichMaze,stats,figstats)

if nargin < 6; sortidx   = [];        end
if nargin < 7; whichMaze = 'accumul'; end

trialType   = cfg.egTrialType;
taxis       = cfg.posBins;
dff         = data.(['ROIavg_' whichMaze '_' trialType]); 
[cl,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.ROIcontraOnly 
  side   = cfg.egTrialType(end);
  idx    = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl = ROIlbl(idx);
  dff    = dff(:,idx,:);
  cl     = cl(idx,:);
end

if cfg.normROImaps
  nt     = size(dff,1);
  minv   = min(dff);
  dff    = (dff - repmat(minv,[nt 1 1]));
  maxv   = max(dff);
  dff    = dff ./ repmat(maxv,[nt 1 1]);
end

dff = mean(dff,3);

if ~isempty(sortidx)
  dff    = dff(:,sortidx);
  ROIlbl = ROIlbl(sortidx);
end

iPanel     = iPanel + 1;
axs        = fig.panel(iPanel);
hold(axs,'on')
for iROI = 1:numel(ROIlbl)
  plot(taxis,dff(:,iROI),'color',cl(iROI,:),'linewidth',.75)
end

if cfg.normPxlMaps
  ylim([0 1])
  set(axs,'ytick',0:.5:1,'xtick',0:50:300)
else
  set(axs,'xtick',0:50:300)
end

xlabel('y pos (cm)')
ylabel('Norm. \DeltaF/F')

% print pvals
if strcmpi(whichMaze,'visGuide')
  text(305,1,'p(task):','color',[.5 .5 .5])
  for iROI = 1:numel(ROIlbl)
    if stats.taskPosANOVA.isSig_task(iROI)
      txtstr = sprintf('%1.2g',stats.taskPosANOVA.pvals_task(iROI));
    else
      txtstr = 'n.s.';
    end
    text(305,.98-.08*iROI,txtstr,'color',cl(iROI,:))
  end
end

switch whichMaze
  case 'visGuide'
    figstats.panelE.dff_norm  = dff;
    figstats.panelE.posBins   = taxis;
    figstats.panelE.ROIlabels = ROIlbl;
    figstats.panelE.taskLabel = 'Vis.-guided';
  otherwise
    figstats.panelF.dff_norm  = dff;
    figstats.panelF.posBins   = taxis;
    figstats.panelF.ROIlabels = ROIlbl;
    figstats.panelE.taskLabel = 'Accum.-towers';
end

end

%% avg ROI activity, lines, zoom in
function iPanel = pnlFramesROIAvgLinesInset(fig,iPanel,data,cfg,ROIlbl,sortidx,whichMaze)

if nargin < 6; sortidx   = [];        end
if nargin < 7; whichMaze = 'accumul'; end

trialType   = cfg.egTrialType;
taxis       = cfg.posBins;
dff         = data.(['ROIavg_' whichMaze '_' trialType]); 

[cl,ROIlbl] = getDefaultROIcl(ROIlbl);

if cfg.ROIcontraOnly 
  side   = cfg.egTrialType(end);
  idx    = cellfun(@(x)(isempty(strfind(x,side))),ROIlbl);
  ROIlbl = ROIlbl(idx);
  dff    = dff(:,idx,:);
  cl     = cl(idx,:);
end

if cfg.normROImaps
  nt     = size(dff,1);
  minv   = min(dff);
  dff    = (dff - repmat(minv,[nt 1 1]));
  maxv   = max(dff);
  dff    = dff ./ repmat(maxv,[nt 1 1]);
end

lastPt      = 150;
idx         = 1:find(taxis==lastPt);
dff         = dff(idx,:,:);
taxis       = taxis(idx);

dff = mean(dff,3);

if ~isempty(sortidx)
  dff    = dff(:,sortidx);
  ROIlbl = ROIlbl(sortidx);
end

iPanel     = iPanel + 1;
axs        = fig.panel(iPanel);
hold(axs,'on')
for iROI = 1:numel(ROIlbl)
  plot(taxis,dff(:,iROI),'color',cl(iROI,:),'linewidth',.75)
end

axis tight

xlabel('y pos (cm)')
ylabel('Norm. \DeltaF/F')

end

%% ANOVA
function stats  = taskANOVA(data,cfg)

whichMaze   = {'accumul','visGuide'};
datavec     = [];
trialType   = cfg.egTrialType;
taxis       = cfg.posBins;

for iMz = 1:numel(whichMaze)
  dff         = data.(['ROIavg_' whichMaze{iMz} '_' trialType]); 
  ROIlbl      = data.ROIlbl;
  [cl,ROIlbl] = getDefaultROIcl(ROIlbl);

  if cfg.ROIcontraOnly 
    side   = cfg.egTrialType(end);
    idx    = cellfun(@(x)(~contains(x,side)),ROIlbl);
    ROIlbl = ROIlbl(idx);
    dff    = dff(:,idx,:);
    cl     = cl(idx,:);
  end

  if cfg.normROImaps
    nt     = size(dff,1);
    minv   = min(dff);
    dff    = (dff - repmat(minv,[nt 1 1]));
    maxv   = max(dff);
    dff    = dff ./ repmat(maxv,[nt 1 1]);
  end
  [nT,nROI,nmice] = size(dff);
  datavec = [datavec; dff(:)];
end

%
roimat   = repmat(ones(nT,nROI)*diag(1:nROI),[1 1 nmice]);
roivec   = repmat(roimat(:),[2 1]);
taskvec  = ones(size(roivec));
taskvec(end/2+1:end) = 2;
mousemat = ones(size(roimat));
for iMouse=1:nmice; mousemat(:,:,iMouse) = mousemat(:,:,iMouse).*iMouse; end
mousevec = repmat(mousemat(:),[2 1]);
timemat  = repmat(repmat(cfg.posBins',[1 nROI]),[1 1 nmice]);
timevec  = repmat(timemat(:),[2 1]);

stats.taskPosANOVA.pvals         = cell(1,nROI);
for iROI = 1:nROI
  stats.taskPosANOVA.pvals{iROI} = anovan(datavec(roivec==iROI), {taskvec(roivec==iROI),timevec(roivec==iROI),mousevec(roivec==iROI)}, ...
                       'model','interaction','varname',{'task','y pos','mouseID'},'display','off');
end
stats.taskPosANOVA.pvals_task = cellfun(@(x)(x(1)),stats.taskPosANOVA.pvals)';
stats.taskPosANOVA.isSig_task = FDR(stats.taskPosANOVA.pvals_task,cfg.alpha);
stats.taskPosANOVA.pvals_pos  = cellfun(@(x)(x(2)),stats.taskPosANOVA.pvals)';
stats.taskPosANOVA.isSig_pos  = FDR(stats.taskPosANOVA.pvals_pos,cfg.alpha);

end