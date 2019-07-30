function Pinto2019_figS4_hemodynamicCorrection(analysisFilePath,summaryFile)

% Pinto2019_figS4_hemodynamicCorrection(analysisFilePath,summaryFile)
% plots quantification of isosbestic correction
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will be saved

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To summarize correction:
%    quantifyVioletCorrection.m (will save BVcorrectionQuantification.mat)(repo: widfieldImaging)
% 2) To generate examples:
%    BVcorrection_paperEgs.m('gp') or ('yfp'). Will save BVcorrection_egFrames.mat (repo: widfieldImaging)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.gcampMice      = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'};
cfg.yfpMice        = {'ty1';'ty2'};
cfg.imageRegion{1} = {101:110,91:100};
cfg.imageRegion{2} = {91:100,76:85};
cfg.egRecGCAMP     = 'ai10/20180409';
cfg.egRecYFP       = 'ty2/20170125';
cfg.egTimePts      = 13001:14200; % 2 min
cfg.frameRate      = 10;
cfg.histBins       = -.1:.005:.1; 
cfg.airpuffEgPath  = 'ai8/20170802_whiskersL_BV/violetCorrection.mat';
cfg.blueCl         = widefieldParams.myblue;
cfg.violetCl       = widefieldParams.mypurple;
cfg.dffCl          = widefieldParams.darkgray;
cfg.yfpCl          = widefieldParams.mediumgreen;
cfg.clLim          = [-.1 .1];

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});

%% load data
data                   = load([analysisFilePath 'BVcorrectionQuantification.mat']);
hemodynCorrect.airpuff = load([analysisFilePath cfg.airpuffEgPath]);

%% compile stats
hemodynCorrect.cfg              = cfg;
hemodynCorrect.stats.recLs      = data.BVrecs;
hemodynCorrect.stats.nGcampMice = numel(data.cfg.gcampMice);
hemodynCorrect.stats.nGcampRecs = sum(cellfun(@(x)(~isempty(x)), ...
                                      reshape(data.gcamp_pxlVals,[1 numel(data.gcamp_pxlVals)])));
hemodynCorrect.stats.nYFPMice   = numel(data.cfg.yfpMice);
hemodynCorrect.stats.nYFPRecs   = sum(cellfun(@(x)(~isempty(x)), ...
                                      reshape(data.yfp_pxlVals,[1 numel(data.yfp_pxlVals)])));
                                    
%% retrieve egs
egData_gcamp = load([analysisFilePath cfg.egRecGCAMP '/BVcorrection_egFrames.mat']);
egData_yfp   = load([analysisFilePath cfg.egRecYFP '/BVcorrection_egFrames.mat']);

hemodynCorrect.taxis                = linspace(0,egData_gcamp.cfg.nFramesTraces/cfg.frameRate-cfg.frameRate,egData_gcamp.cfg.nFramesTraces);

hemodynCorrect.gcamp_eg.dff_blueExc = egData_gcamp.dffReg.blue;
hemodynCorrect.gcamp_eg.dff_violExc = egData_gcamp.dffReg.violet; 
hemodynCorrect.gcamp_eg.dff_correct = egData_gcamp.dffReg.correct;
hemodynCorrect.gcamp_eg.cfg         = egData_gcamp.cfg;
hemodynCorrect.gcamp_eg.frames      = egData_gcamp.frames;

hemodynCorrect.yfp_eg.dff_blueExc   = egData_yfp.dffReg.blue;
hemodynCorrect.yfp_eg.dff_violExc   = egData_yfp.dffReg.violet; 
hemodynCorrect.yfp_eg.dff_correct   = egData_yfp.dffReg.correct;
hemodynCorrect.yfp_eg.cfg           = egData_yfp.cfg;
hemodynCorrect.yfp_eg.frames        = egData_yfp.frames;

%% build histograms 
for iReg = 1:numel(cfg.imageRegion)
  gcamp  = []; 
  gcampb = [];
  gcampv = [];
  yfp    = [];
  yfpb   = [];
  yfpv   = [];
  for iMouse = 1:numel(data.cfg.gcampMice)
    for iRec = 1:data.cfg.recsPerAnimal
      gcamp  = [gcamp; data.gcamp_pxlVals{iMouse,iRec,iReg}];
      gcampb = [gcampb; data.gcamp_pxlValsB{iMouse,iRec,iReg}];
      gcampv = [gcampv; data.gcamp_pxlValsV{iMouse,iRec,iReg}];
    end
  end
  for iMouse = 1:numel(data.cfg.yfpMice)
    for iRec = 1:data.cfg.recsPerAnimal
      yfp    = [yfp; data.yfp_pxlVals{iMouse,iRec,iReg}];
      yfpb   = [yfpb; data.yfp_pxlValsB{iMouse,iRec,iReg}];
      yfpv   = [yfpv; data.yfp_pxlValsV{iMouse,iRec,iReg}];
    end
  end
  hemodynCorrect.hist_xaxis{iReg}      = toBinCenters(cfg.histBins);
  hemodynCorrect.hist_gcamp{iReg}      = histcounts(gcamp,cfg.histBins,'Normalization','probability');
  hemodynCorrect.hist_gcamp_blue{iReg} = histcounts(gcampb,cfg.histBins,'Normalization','probability');
  hemodynCorrect.hist_gcamp_viol{iReg} = histcounts(gcampv,cfg.histBins,'Normalization','probability');
  hemodynCorrect.hist_yfp{iReg}        = histcounts(yfp,cfg.histBins,'Normalization','probability');
  hemodynCorrect.hist_yfp_blue{iReg}   = histcounts(yfpb,cfg.histBins,'Normalization','probability');
  hemodynCorrect.hist_yfp_viol{iReg}   = histcounts(yfpv,cfg.histBins,'Normalization','probability');
end

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'hemodynCorrect','-v7.3')
  else
    save(summaryFile,'hemodynCorrect','-append')
  end
end

%% Configure figure and panels
layout       = [ 1  2  3  4  5  6  7  8  ...
               ; 9  10 11 12 13 14 15 16 ...
               ; 17 17 17 17 19 19 19 19 ...
               ; 18 18 18 18 20 20 20 20 ...
               ; 21 22 23 24 25 26 27 28 ...
               ; 29 30 31 32 33 34 35 36 ...
               ; 37 37 37 37 39 39 39 39 ...
               ; 38 38 38 38 40 40 40 40 ...
               ; 41 42 43 44 45 46 47 48 ...
               ];
fig          = PaneledFigure(layout, 'tiny');
iPanel       = 0;
stats        = [];

%% Plot frames gcamp 
iPanel       = pnlFrames(fig, iPanel, hemodynCorrect, cfg, 'gcamp');

%% Plot blue and violet, dff gcamp (V1)
iPanel       = pnlTrace(fig, iPanel, hemodynCorrect, cfg, 'gcamp', 1);

%% Plot blue and violet, dff gcamp (Retrosplenial)
iPanel       = pnlTrace(fig, iPanel, hemodynCorrect, cfg, 'gcamp', 2);

%% Plot frames yfp 
iPanel       = pnlFrames(fig, iPanel, hemodynCorrect, cfg, 'yfp');

%% Plot blue and violet, dff yfp (V1)
iPanel       = pnlTrace(fig, iPanel, hemodynCorrect, cfg, 'yfp', 1);

%% Plot blue and violet, dff gcamp (Retrosplenial)
iPanel       = pnlTrace(fig, iPanel, hemodynCorrect, cfg, 'yfp', 2);

%% Plot histogram
iPanel       = pnlHist(fig, iPanel, hemodynCorrect, cfg, 'gcamp', 1);

%% Plot histogram
iPanel       = pnlHist(fig, iPanel, hemodynCorrect, cfg, 'yfp', 1);

%% Plot airpuff
iPanel       = pnlAirpuff(fig, iPanel, hemodynCorrect.airpuff, cfg);

%% Save figure to disk
versionInfo.stats = hemodynCorrect.stats;
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% blue /dff eg frames
function iPanel  = pnlFrames(fig, iPanel, data, cfg, whichData)

% frames
B         = data.([whichData '_eg']).frames.dffBlue;
C         = data.([whichData '_eg']).frames.dff;
cmap      = red2blue;
cmap(1,:) = [0 0 0];
[nx,~,nz] = size(B);
barsize   = widefieldParams.pxlPerMM / widefieldParams.dsFactor;

% ROI squares
ROI1      = data.([whichData '_eg']).cfg.imageRegion{1};
ROI2      = data.([whichData '_eg']).cfg.imageRegion{2};
xvec1     = [ROI1{2}(1) ROI1{2}(end) ROI1{2}(end) ROI1{2}(1) ROI1{2}(1)];
yvec1     = nx - [ROI1{1}(1) ROI1{1}(1) ROI1{1}(end) ROI1{1}(end) ROI1{1}(1)];
xvec2     = [ROI2{2}(1) ROI2{2}(end) ROI2{2}(end) ROI2{2}(1) ROI2{2}(1)];
yvec2     = nx - [ROI2{1}(1) ROI2{1}(1) ROI2{1}(end) ROI2{1}(end) ROI2{1}(1)];

% blue excit.
for iFrame = 1:nz
  iPanel     = iPanel + 1;
  axs        = fig.panel(iPanel);
  hold(axs,'on')
  imagesc(flipud(B(:,:,iFrame)),cfg.clLim)
  axis image; axis off
  colormap(cmap)
  
  if iFrame == 1; plot([5 5+barsize],[5 5],'w-','linewidth',2.5); end% scale bar
  
  % squares
  h = fill(xvec1,yvec1,[.8 .8 .8]);
  set(h,'facecolor','none','edgecolor',[.8 .8 .8],'linewidth',.75);
  h = fill(xvec2,yvec2,[.8 .8 .8]);
  set(h,'facecolor','none','edgecolor',[.8 .8 .8],'linewidth',.75);
  
  if iFrame == round(nz/2)
    title([whichData ', 470-nm excit.'])
  end
end
cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\DeltaF/F';

% corrected
for iFrame = 1:nz
  iPanel     = iPanel + 1;
  axs        = fig.panel(iPanel);
  hold(axs,'on')
  imagesc(flipud(C(:,:,iFrame)),cfg.clLim)
  axis image; axis off
  colormap(cmap)
  
  if iFrame == 1; plot([5 5+barsize],[5 5],'w-','linewidth',2.5); end% scale bar
  
  % squares
  h = fill(xvec1,yvec1,[.8 .8 .8]);
  set(h,'facecolor','none','linewidth',.75);
  h = fill(xvec2,yvec2,'k');
  set(h,'facecolor','none','linewidth',.75);
  
  if iFrame == round(nz/2)
    title([whichData ', corrected'])
  end
end
cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = '\DeltaF/F';

end

%% ------------------------------------------------------------------------
%% blue / violet / dff eg
function iPanel  = pnlTrace(fig, iPanel, data, cfg, whichData, whichRegion)

B = data.([whichData '_eg']).dff_blueExc{whichRegion};
V = data.([whichData '_eg']).dff_violExc{whichRegion};
C = data.([whichData '_eg']).dff_correct{whichRegion};

switch whichRegion
  case 1
    regLbl = 'Vis. ctx';
  case 2
    regLbl = 'RSC';
end

% blue violet
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

% shading for frame eg period
nfim = data.([whichData '_eg']).cfg.nFramesIm*data.([whichData '_eg']).cfg.frameStep;
xvec = [data.taxis(floor(end/2)-nfim/2) data.taxis(floor(end/2)+nfim/2) ...
        data.taxis(floor(end/2)+nfim/2) data.taxis(floor(end/2)-nfim/2) data.taxis(floor(end/2)-nfim/2)];
yvec = [cfg.clLim(1) cfg.clLim(1) cfg.clLim(2) cfg.clLim(2) cfg.clLim(1)];
h = fill(xvec,yvec,[.8 .8 .8]);
set(h,'edgecolor','none','facealpha',.5);

% data
plot([data.taxis(1) data.taxis(end)],[0 0],'--','color',FormatDefaults.lightGray)
h1 = plot(data.taxis,B,'-','color',cfg.blueCl,'linewidth',1);
h2 = plot(data.taxis,V,'-','color',cfg.violetCl,'linewidth',1);

legend([h1(1); h2(1)],{'470-nm excit.';'410-nm excit.'},'location','best'); legend('boxoff')
ylim(cfg.clLim); 
set(axs,'xtick',[],'xcolor','w','ytick',-.1:.1:.1)
ylabel('\DeltaF/F')
title([whichData ' - ' regLbl],'fontsize',13)

% corrected
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

% shading for frame eg period
h = fill(xvec,yvec,[.8 .8 .8]);
set(h,'edgecolor','none','facealpha',.5);

plot([data.taxis(1) data.taxis(end)],[0 0],'--','color',FormatDefaults.lightGray)
h1 = plot(data.taxis,C,'-','color',cfg.dffCl,'linewidth',1);

legend(h1(1),{'corrected'},'location','best'); legend('boxoff')
ylim(cfg.clLim); 
set(axs,'xtick',[],'xcolor','w','ytick',-.1:.1:.1)
ylabel('\DeltaF/F')
% title('GCaMP6f','fontsize',13)

% time scale
plot([0 10],[cfg.clLim(1) cfg.clLim(1)],'k-')
text(5, -.082, '10 s', 'horizontalAlignment', 'center', 'color', 'k', 'fontsize', 10)

end

%% ------------------------------------------------------------------------
%% histogram dff values
function iPanel  = pnlHist(fig, iPanel, data, cfg, whichData, whichRegion)

xaxis = data.hist_xaxis{whichRegion};
B     = data.(['hist_' whichData '_blue']){whichRegion};
V     = data.(['hist_' whichData '_viol']){whichRegion};
C     = data.(['hist_' whichData ]){whichRegion};

% blue excitation only
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

hHisto                  = gobjects(0);
hHisto(end+1)           = bar ( axs, xaxis, B, FormatDefaults.barWidth                 ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , widefieldParams.myblue                 ...
                              );
hHisto(end+1)           = bar ( axs, xaxis, V, FormatDefaults.barWidth             ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , widefieldParams.mypurple                 ...
                              );

% Simulate transparency (actual transparency has PDF output issues)
colors                  = [widefieldParams.myblue; widefieldParams.mypurple];
                          bar ( axs, xaxis, min([B; V],[],1), FormatDefaults.barWidth ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , mean(colors,1)                            ...
                              );

xlabel(axs, '\DeltaF/F'); 
ylabel(axs, 'Prop. frames'); 
set(axs,'yscale','log')
% axis tight
xlim([cfg.histBins(1)-mode(diff(cfg.histBins)) cfg.histBins(end)+mode(diff(cfg.histBins))])
ylim([1e-4 .5])      
set(axs,'ytick',[1e-4 1e-2])

legend( hHisto, {'470 nm','410 nm'}                               ...
      , 'Location'      , 'NorthEast'                             ...
      , 'FontSize'      , FormatDefaults.legendFontSize           ...
      , 'Box'           , 'off'                                   ...
      );
      
    
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs,'on')

hHisto                  = gobjects(0);
hHisto(end+1)           = bar ( axs, xaxis, C, FormatDefaults.barWidth                 ...
                              , 'EdgeColor'     , 'none'                                    ...
                              , 'FaceColor'     , FormatDefaults.lightGray                 ...
                              );

legend( hHisto, {'corrected'}                               ...
      , 'Location'      , 'NorthEast'                             ...
      , 'FontSize'      , FormatDefaults.legendFontSize           ...
      , 'Box'           , 'off'                                   ...
      );
      
    
xlabel(axs, '\DeltaF/F'); 
ylabel(axs, 'Prop. frames'); 
set(axs,'yscale','log')
axis tight
xlim([cfg.histBins(1)-mode(diff(cfg.histBins)) cfg.histBins(end)+mode(diff(cfg.histBins))])
ylim([1e-4 .5])      
set(axs,'ytick',[1e-4 1e-2])

end

%% ------------------------------------------------------------------------
%% correction airpuff expt
function iPanel  = pnlAirpuff(fig, iPanel, data, cfg)

iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);

hold(axs,'on')
xl = [-1 5];
yl = [-.02 .02];
plot(xl,[0 0],'--','color',FormatDefaults.lightGray);
fh = fill([0 .2 .2 0],[yl(1) yl(1) yl(2) yl(2)],FormatDefaults.lightGray);
fh.EdgeColor = FormatDefaults.lightGray;
fh.FaceAlpha = .5;
h1 = plot(data.taxis, data.dffMeanV, '-', 'color', cfg.violetCl, 'linewidth', 1);
h2 = plot(data.taxis, data.dffMeanB, '-', 'color', cfg.blueCl,   'linewidth', 1);
h3 = plot(data.taxis, data.dffMean,  '-', 'color', cfg.dffCl,    'linewidth', 1);
plot(data.taxis, data.dffMeanV - data.dffSemV, '--', 'color', cfg.violetCl, 'linewidth', .5)
plot(data.taxis, data.dffMeanB - data.dffSemB, '--', 'color', cfg.blueCl,   'linewidth', .5)
plot(data.taxis, data.dffMean  - data.dffSem,  '--', 'color', cfg.dffCl,    'linewidth', .5)
plot(data.taxis, data.dffMeanV + data.dffSemV, '--', 'color', cfg.violetCl, 'linewidth', .5)
plot(data.taxis, data.dffMeanB + data.dffSemB, '--', 'color', cfg.blueCl,   'linewidth', .5)
plot(data.taxis, data.dffMean  + data.dffSem,  '--', 'color', cfg.dffCl,    'linewidth', .5)
legend([h1(1); h2(1); h3(1)],{'405','473','corrected'},'location','best'); legend('boxoff')
axis tight; 
set(axs,'ytick',-.02:.01:.02)
xlabel('Time from airpuff (s)')
ylabel('\DeltaF/F')

end