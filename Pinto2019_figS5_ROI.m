function Pinto2019_figS5_ROI(analysisFilePath,summaryFile)

% Pinto2019_figS5_ROI(analysisFilePath,summaryFile)
% plots quantification of isosbestic correction
% analysisFilePath is path for data analysis files to be loaded
% summaryFile is path for file where summary stats will saved
% this function does not require pre-saved summary data, just ROI and
% sensory maps previously generated as described in Pinto2019_fig3_WFdynamics.m

%% Analysis configuration
cfg.egMap          = 'ai3';
cfg.egROI          = {'VISp','mV2','VISa','SS'};
cfg.egROIlbl       = {'V1','mV2','PPC','SS'};
cfg.ROIcl          = widefieldParams.areaCl;
cfg.ROIEgCl        = cfg.ROIcl([1 2 3 8],:);
cfg.ROIorigName    = {'VISp','mV2','VISa','RSP','aMOs','mMOs','MOp','SS'}; 
cfg.ROIlbl         = {'V1','mV2','PPC','RSC','aM2','mM2','M1','SS'};
cfg.ROIlblOrder    = [1:4 6 5 7 8];
cfg.ROIorder       = [1:8 11 12 9 10 13:16];

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% compile maps
regMaps            = registerMapsToAllen([analysisFilePath 'maps_all.mat']);

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'regMaps','-v7.3')
  else
    save(summaryFile,'regMaps','-append')
  end
end
%% Configure figure and panels
layout       = [  1  3  4  5  6   ...
                ; 2  7  8  9  10  ...
               ];
fig          = PaneledFigure(layout, 'small');
iPanel       = 0;

%% Plot colored ROI outlines
iPanel       = pnlEgROI(fig, iPanel, regMaps, cfg);

%% Plot example GCAMP map (visual, airpuff w/ V1, PPC, SS outlines)
plotEg       = true;
iPanel       = pnlMap(fig, iPanel, cfg, regMaps, plotEg);

%% Plot avg GCAMP map (visual, airpuff w/ V1, PPC, SS outlines)
plotEg       = false;
iPanel       = pnlMap(fig, iPanel, cfg, regMaps, plotEg);

%% Save figure to disk
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])

end

%% ------------------------------------------------------------------------
%% plot ROI on eg brain 
function iPanel = pnlEgROI(fig,iPanel,regMaps,cfg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

%%
idx    = strcmpi(regMaps.gcamp.mice,cfg.egMap);
im     = uint8(regMaps.gcamp.meanproj{idx});
rois   = regMaps.gcamp.outline{idx};
rois   = rois(cfg.ROIorder);

image(im); axis off; axis image
ct     = 1; hold on
for iROI = 1:2:numel(rois)
  plot(rois{iROI}(:,2),rois{iROI}(:,1),'linewidth',2,'color',cfg.ROIcl(ct,:));
  plot(rois{iROI+1}(:,2),rois{iROI+1}(:,1),'linewidth',2,'color',cfg.ROIcl(ct,:));
  ct = ct+1;
end
set(axs,'view',[180 90]);

% scale bar
ds    = size(im,1)/widefieldParams.imsize(1);
scale = widefieldParams.pxlPerMM*ds;
plot([10 10+scale],[120 120],'w-','linewidth',2)
text(10+scale/2,110,'1 mm','fontsize',10,'color','w','horizontalAlignment','center')
%%
% ROI lbls
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')
xlim([0 1]); ylim([0 1]); axis off

lbl = cfg.ROIlbl(cfg.ROIlblOrder);
for iROI = 1:numel(cfg.ROIlbl)
  text(.1,iROI*.12,lbl{iROI},'color',cfg.ROIcl(iROI,:),'fontsize',14)
end

end

%% ------------------------------------------------------------------------
%% plot ROI on eg brain 
function iPanel = pnlMap(fig,iPanel,cfg,regMaps,plotEg)

iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')

%% map
if plotEg
  idx       = strcmpi(regMaps.gcamp.mice,cfg.egMap);
  fieldSign = regMaps.gcamp.fieldSign{idx};
  airpuff   = regMaps.gcamp.airpuff{idx};
  rois      = regMaps.gcamp.outline{idx}; 
else
  fieldSign = regMaps.gcamp_avg.fieldSign;
  airpuff   = regMaps.gcamp_avg.airpuff;
  rois      = regMaps.gcamp_avg.outline; 
end

%% retrieve only relevant ROIs
egrois = {};
roiid  = [];
for iEg = 1:numel(cfg.egROI)
  isROI = find(cellfun(@(x)(~isempty(strfind(x,cfg.egROI{iEg}))),regMaps.gcamp_avg.ROIlbl));
  egrois(end+1:end+numel(isROI)) = rois(isROI);
  roiid(end+1:end+numel(isROI))  = iEg;
end

%% make combined RGB image
im1 = ones(size(airpuff,1),size(airpuff,2));
im2 = im1;
im3 = im1;

% red scale (negative field sign)
im1(fieldSign<0) = 1;
im2(fieldSign<0) = 1+fieldSign(fieldSign<0);
im3(fieldSign<0) = 1+fieldSign(fieldSign<0);

% blue scale (positive field sign)
im3(fieldSign>0) = 1;
im2(fieldSign>0) = 1-fieldSign(fieldSign>0);
im1(fieldSign>0) = 1-fieldSign(fieldSign>0);

im               = cat(3,cat(3,im1,im2),im3);
im               = uint8(im.*255);

%% plot
image(fliplr(im)); axis off; axis image; hold on
for iROI = 1:numel(cfg.egROI)
  idx = find(roiid == iROI);
  for iSide = 1:numel(idx)
    plot(egrois{idx(iSide)}(:,2),egrois{idx(iSide)}(:,1),'linewidth',2,'color',cfg.ROIEgCl(iROI,:));
  end
end
set(axs,'view',[180 90]);

%% color bars
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')
xlim([0 1]); ylim([0 1]); axis off

if plotEg
  text(1,1,'Example mouse','color',[.5 .5 .5],'fontsize',12,'horizontalAlignment','center')
  red   = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)'];
  blue  = fliplr(red);

  colors     = [red;flipud(blue)];
  colors     = reshape(colors(:,:),[],1,3);
  image ( 'Parent'              , axs              ...
        , 'XData'               , [.1 .125]        ...
        , 'YData'               , [.1 .5]          ...
        , 'CData'               , colors           ...
        );
  t1 = text(.2,.1,'-1','fontsize',10,'color','k','horizontalAlignment','center');
  set(t1,'rotation',90);
  t2 = text(.2,.3,'0','fontsize',10,'color','k','horizontalAlignment','center');
  set(t2,'rotation',90);
  t3 = text(.2,.5,'1','fontsize',10,'color','k','horizontalAlignment','center');
  set(t3,'rotation',90);
  t4 = text(.3,.3,'Visual field sign','fontsize',12,'color','k','horizontalAlignment','center');
  set(t4,'rotation',90);
else
  text(1,1,sprintf('Average (n = %d)',numel(regMaps.gcamp.mice)), ...
    'color',[.5 .5 .5],'fontsize',12,'horizontalAlignment','center')
end

%% next panel airpuff
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')
  
% magenta scale (airpuff response)
im1        = ones(size(airpuff,1),size(airpuff,2));
im2        = im1;
im3        = im1;
ima        = airpuff ./ max(airpuff(:));
im2(ima>0) = 1-ima(ima>0);
% im3(ima>0) = 1-ima(ima>0);
% im1(ima>0) = 1-ima(ima>0);
im         = cat(3,cat(3,im1,im2),im3);
im         = uint8(im.*255);

%% plot
image(fliplr(im)); axis off; axis image; hold on
for iROI = 1:numel(cfg.egROI)
  idx = find(roiid == iROI);
  for iSide = 1:numel(idx)
    plot(egrois{idx(iSide)}(:,2),egrois{idx(iSide)}(:,1),'linewidth',2,'color',cfg.ROIEgCl(iROI,:));
  end
end
set(axs,'view',[180 90]);


iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs,'on')
xlim([0 1]); ylim([0 1]); axis off

if plotEg  
  magenta    = flipud([ones(100,1) linspace(0,1,100)' ones(100,1)]);
  colors     = reshape(magenta(:,:),[],1,3);
  image ( 'Parent'              , axs              ...
        , 'XData'               , [.1 .125]        ...
        , 'YData'               , [.1 .5]          ...
        , 'CData'               , colors           ...
        );
  t1 = text(.2,.1,'0','fontsize',10,'color','k','horizontalAlignment','center');
  set(t1,'rotation',90);
  t3 = text(.2,.5,'1','fontsize',10,'color','k','horizontalAlignment','center');
  set(t3,'rotation',90);
  t4 = text(.3,.3,sprintf('Face pad airpuff resp.\n(Norm. \DeltaF/F)'),'fontsize',12,'color','k','horizontalAlignment','center');
  set(t4,'rotation',90);

end
end