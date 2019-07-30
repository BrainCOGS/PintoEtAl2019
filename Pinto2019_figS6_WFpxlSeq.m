function Pinto2019_figS6_WFpxlSeq(analysisFilePath,summaryFile)

% Pinto2019_figS6_WFpxlSeq(analysisFilePath,summaryFile)
% plots examples of pxl activation sequences
% analysisFilePath is path for data analysis files to be loaded

%% ------------------------------------------------------------------------
%% Instructions for running analysis from raw data
% 1) To summarize experiments and generate sorted pixel matrices:
%    summarizeAvgDff.m (will save avgDffSummary.mat, loaded here)(repo: widefieldImaging)
%    (calls activitySequences.m)
%% ------------------------------------------------------------------------

%% Analysis configuration
cfg.egRec            = 'ai2/20170125';
cfg.compressPctile   = [1 99]; % compress dff scale in pxl maps for display
cfg.seqType          = 'pxl'; % show quantification of pxl or ROI COM sequences
cfg.corrType         = 'COM'; % OVERALL or COM

%% code version
fprintf('collecting version info...\n')
codeFile           = mfilename('fullpath');
versionInfo        = collectVersionInfo(codeFile, cfg, [], [], {});
cfg.stockFigPath   = [getRepositoryPath(codeFile) '/stockFigs'];

%% load data
fprintf('loading summary data from disk...\n')
load([analysisFilePath '/avgDffSummary.mat'],'avgDffSumm');

%% get desired example
mouse = cfg.egRec(1:3);
rec   = cfg.egRec(5:end);
midx  = strcmp(avgDffSumm.cfg.mice,mouse);
ridx  = cellfun(@(x)(~isempty(strfind(x,rec))),avgDffSumm.mouse(midx).recls);
seq   = avgDffSumm.mouse(midx).pxl.seq.rec(ridx);

%% get pixel identities for ROI color coding
load([analysisFilePath '/' cfg.egRec '/dff.mat'],'dff');
load([analysisFilePath '/' cfg.egRec '/ROIfromRef.mat'],'ROIlbl','ROI');
mask            = dff(:,:,1); clear dff
nanMask         = isnan(mask);
[nX,nY]         = size(nanMask);
rowID           = repmat((1:nX)',[1 nY]);
colID           = repmat(1:nY,[nX 1]);
rowID(nanMask)  = nan;
colID(nanMask)  = nan;
rowID           = conditionDffMat(cat(3,rowID,rowID));
rowID           = rowID(1,:);
colID           = conditionDffMat(cat(3,colID,colID));
colID           = colID(1,:);
nROI            = numel(rowID);

%% Configure figure and panels
layout       = [ 1:4  ...
               ; 5:8  ...
               ; 9:12 ...
               ];
fig          = PaneledFigure(layout, 'tall');
iPanel       = 0;

%% maze schematic
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematic.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image

%% towers sequences 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

idx  = seq.accumul.R.sortIdx;
mat  = seq.accumul.R.sortedbyR.dff;
mat  = mat - repmat(min(mat),[size(mat,1) 1]);
mat  = mat ./ repmat(max(mat),[size(mat,1) 1]);
imagesc(seq.bins,1:nROI,mat',[0 1])

% plot colored lines indicating pxl ROI membership
roiIm = zeros(nROI,40,3);
for iPxl = 1:nROI
  cl = getPxlROIcl([rowID(idx(iPxl)) colID(idx(iPxl))],ROI,ROIlbl);
  roiIm(iPxl,:,:) = repmat(reshape(cl,[1 1 3]),[1 40 1]);
end
xlim([0 seq.bins(end)+60])
ylim([1 nROI])
image(axs,[seq.bins(end)+21 seq.bins(end)+60],[1 nROI],roiIm)
set(axs, 'ytick', []);
axis ij
xlabel('y (cm)'); ylabel('Pixel'); 
text(0,-200,'towers, R choice, sorted by R','color',[.6 .6 .6],'fontsize',FormatDefaults.legendFontSize)

%%
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

idx  = seq.accumul.L.sortIdx;
mat  = seq.accumul.R.sortedbyL.dff;
mat  = mat - repmat(min(mat),[size(mat,1) 1]);
mat  = mat ./ repmat(max(mat),[size(mat,1) 1]);
imagesc(seq.bins,1:nROI,mat',[0 1])

% plot colored lines indicating pxl ROI membership
roiIm = zeros(nROI,40,3);
for iPxl = 1:nROI
  cl = getPxlROIcl([rowID(idx(iPxl)) colID(idx(iPxl))],ROI,ROIlbl);
  roiIm(iPxl,:,:) = repmat(reshape(cl,[1 1 3]),[1 40 1]);
end
xlim([0 seq.bins(end)+60])
ylim([1 nROI])
image(axs,[seq.bins(end)+21 seq.bins(end)+60],[1 nROI],roiIm)
set(axs, 'ytick', []);
axis ij
xlabel('y (cm)'); ylabel('Pixel'); 
text(0,-200,'towers, R choice, sorted by L','color',[.6 .6 .6],'fontsize',FormatDefaults.legendFontSize)

%% captions 
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
xlim([0 1]); ylim([0 1])
text(axs,.05,.9,sprintf('COM corr (accumul, R vs L) = %1.2f\nCOM corr (visGuide, R vs L) = %1.2f\nCOM corr (accumul vs visGuide) = %1.2f',...
  seq.accumul.COMcorr_RvsL,seq.visGuide.COMcorr_RvsL,seq.COMcorr_accumulVSvisGuide),...
  'color','k','fontsize',11)

lbls = unique(cellfun(@(x)(x(1:end-2)),ROIlbl,'UniformOutput',false));
for iROI = 1:numel(lbls)
  [cl,lbl] = getDefaultROIcl(lbls{iROI});
  text(.05,.7-(iROI-1)*.06,lbl,'color',cl,'fontsize',11)
end
text(axs,.05,.7-iROI*.06,'Not assigned to any ROI','color',[.5 .5 .5],'fontsize',11)

%% maze schematic
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
im     = imread([cfg.stockFigPath '/mazeSchematicVisualGuide.png']);
[nX,~] = size(im);

hold(axs, 'on')
image(flipud(im)); axis off; axis image

%% visual guide sequences
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

idx  = seq.visGuide.R.sortIdx;
mat  = seq.visGuide.R.sortedbyR.dff;
mat  = mat - repmat(min(mat),[size(mat,1) 1]);
mat  = mat ./ repmat(max(mat),[size(mat,1) 1]);
imagesc(seq.bins,1:nROI,mat',[0 1])

% plot colored lines indicating pxl ROI membership
roiIm = zeros(nROI,40,3);
for iPxl = 1:nROI
  cl = getPxlROIcl([rowID(idx(iPxl)) colID(idx(iPxl))],ROI,ROIlbl);
  roiIm(iPxl,:,:) = repmat(reshape(cl,[1 1 3]),[1 40 1]);
end
xlim([0 seq.bins(end)+60])
ylim([1 nROI])
image(axs,[seq.bins(end)+21 seq.bins(end)+60],[1 nROI],roiIm)
set(axs, 'ytick', []);
axis ij
xlabel('y (cm)'); ylabel('Pixel'); 
text(0,-200,'warm-up, R choice, sorted by R','color',[.6 .6 .6],'fontsize',FormatDefaults.legendFontSize)

%%
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

idx  = seq.visGuide.L.sortIdx;
mat  = seq.visGuide.R.sortedbyL.dff;
mat  = mat - repmat(min(mat),[size(mat,1) 1]);
mat  = mat ./ repmat(max(mat),[size(mat,1) 1]);
imagesc(seq.bins,1:nROI,mat',[0 1])

% plot colored lines indicating pxl ROI membership
roiIm = zeros(nROI,40,3);
for iPxl = 1:nROI
  cl = getPxlROIcl([rowID(idx(iPxl)) colID(idx(iPxl))],ROI,ROIlbl);
  roiIm(iPxl,:,:) = repmat(reshape(cl,[1 1 3]),[1 40 1]);
end
xlim([0 seq.bins(end)+60])
ylim([1 nROI])
image(axs,[seq.bins(end)+21 seq.bins(end)+60],[1 nROI],roiIm)
set(axs, 'ytick', []);
axis ij
xlabel('y (cm)'); ylabel('Pixel'); 
text(0,-200,'warmu-up, R choice, sorted by L','color',[.6 .6 .6],'fontsize',FormatDefaults.legendFontSize)


%%
iPanel = iPanel + 1;
axs    = fig.panel(iPanel);
hold(axs, 'on')

idx  = seq.accumul.R.sortIdx;
mat  = seq.visGuide.R.sortedbyRaccumul.dff;
mat  = mat - repmat(min(mat),[size(mat,1) 1]);
mat  = mat ./ repmat(max(mat),[size(mat,1) 1]);
imagesc(seq.bins,1:nROI,mat',[0 1])

% plot colored lines indicating pxl ROI membership
roiIm = zeros(nROI,40,3);
for iPxl = 1:nROI
  cl = getPxlROIcl([rowID(idx(iPxl)) colID(idx(iPxl))],ROI,ROIlbl);
  roiIm(iPxl,:,:) = repmat(reshape(cl,[1 1 3]),[1 40 1]);
end
xlim([0 seq.bins(end)+60])
ylim([1 nROI])
image(axs,[seq.bins(end)+21 seq.bins(end)+60],[1 nROI],roiIm)
set(axs, 'ytick', []);
axis ij
xlabel('y (cm)'); ylabel('Pixel'); 
text(0,-200,'warmp-up, R choice, sorted by towers R','color',[.6 .6 .6],'fontsize',FormatDefaults.legendFontSize);

cbar              = smallcolorbar(axs);
pos               = cbar.Position;
cbar.Units        = get(axs,'units');
cbar.Position     = pos;
cbar.Label.String = 'Norm. \DeltaF/F';

[iPanel,pxlSeqStats]  = pnlSeqSimilarity(fig, iPanel, avgDffSumm, cfg);

%% save analysis summary file
if ~isempty(summaryFile)
  if isempty(dir(summaryFile))
    save(summaryFile,'pxlSeqStats','-v7.3')
  else
    save(summaryFile,'pxlSeqStats','-append')
  end
end

%% Save figure to disk
fig.export(codeFile, versionInfo, true, true);
delete([codeFile '/panel*'])
end

%% find pxl ROI membership and return color
function cl = getPxlROIcl(coord,ROI,ROIlbl)

isROI = cellfun(@(x)(sum(x(:,1)==coord(1) & x(:,2)==coord(2))>0),ROI);
if sum(isROI) == 0
  cl = [.5 .5 .5];
else
	cl = getDefaultROIcl(ROIlbl(isROI)); 
end

end

%% ------------------------------------------------------------------------
%% bar plots comparing similarity of COMs for different conditions
function [iPanel,stats] = pnlSeqSimilarity(fig, iPanel, avgDffSumm, cfg, stats)

iPanel           = iPanel + 1;
axs              = fig.panel(iPanel);
hold(axs,'on')

thismean = mean(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumul_RvsL']));
ndata    = numel(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumul_RvsL']));
thissem  = std(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumul_RvsL'])) ./ sqrt(ndata-1);
bar(1,thismean,'edgecolor',widefieldParams.mediumgray,'facecolor',widefieldParams.mediumgray)
errorbar(1,thismean,thissem,'-','color',widefieldParams.mediumgray)

thismean = mean(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_visGuide_RvsL']));
thissem  = std(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_visGuide_RvsL'])) ./ sqrt(ndata-1);
bar(2,thismean,'edgecolor',widefieldParams.mediumgray,'facecolor',widefieldParams.mediumgray)
errorbar(2,thismean,thissem,'-','color',widefieldParams.mediumgray)

thismean = mean(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumulVSvisGuide']));
thissem  = std(avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumulVSvisGuide'])) ./ sqrt(ndata-1);
bar(3,thismean,'edgecolor',widefieldParams.mediumgray,'facecolor',widefieldParams.mediumgray)
errorbar(3,thismean,thissem,'-','color',widefieldParams.mediumgray)

set(gca,'xtick',1:3,'xticklabel',{'towers task, R vs L','ctrl task, R vs L','towers vs ctrl'})
rotateXLabels(gca,75);
xlim([.25 3.75]);
ylabel('C.O.M. corr. (r)')

%% ANOVA and post hoc
vals = [avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumul_RvsL']);      ...
        avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_visGuide_RvsL']);     ...
        avgDffSumm.([cfg.seqType 'Seq_' cfg.corrType 'corr_accumulVSvisGuide']); ...
        ];
mid  = repmat((1:ndata)',[3 1]);
sqid = [ones(ndata,1); ones(ndata,1)+1; ones(ndata,1)+2];

[stats.COMcorrANOVA.p,stats.COMcorrANOVA.table,stats.COMcorrANOVA.stats] ...
               = anovan(vals,{sqid,mid},'varnames',{'type';'mice'},'display','off');
stats.multcomp = multcompare(stats.COMcorrANOVA.stats,'display','off');
p_1vs2         = stats.multcomp(1,end);
p_1vs3         = stats.multcomp(2,end);
p_2vs3         = stats.multcomp(3,end);

% plot
yl = get(axs,'ylim');
plot([1.1 2.9],[yl(2)*.99 yl(2)*.99],'-','color',[.3 .3 .3])
plot([1.1 1.9],[yl(2)*.95 yl(2)*.95],'-','color',[.3 .3 .3])
plot([1.9 2.9],[yl(2)*.83 yl(2)*.83],'-','color',[.3 .3 .3])

if p_1vs3 >= .05; t1 = 'n.s.'; end
if p_1vs3 < .05;  t1 = '*';    end
if p_1vs3 < .01;  t1 = '**';   end
if p_1vs3 < .001; t1 = '***';  end
text(2,yl(2)*1.05,t1,'fontsize',8,'color',[.3 .3 .3])

if p_1vs2 >= .05; t1 = 'n.s.'; end
if p_1vs2 < .05;  t1 = '*';    end
if p_1vs2 < .01;  t1 = '**';   end
if p_1vs2 < .001; t1 = '***';  end
text(1.5,yl(2),t1,'fontsize',8,'color',[.3 .3 .3])

if p_2vs3 >= .05; t1 = 'n.s.'; end
if p_2vs3 < .05;  t1 = '*';    end
if p_2vs3 < .01;  t1 = '**';   end
if p_2vs3 < .001; t1 = '***';  end
text(2.5,yl(2)*.94,t1,'fontsize',8,'color',[.3 .3 .3])

end
