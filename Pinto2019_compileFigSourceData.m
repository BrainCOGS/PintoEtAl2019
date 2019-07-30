function Pinto2019_compileFigSourceData(figDataFile,summaryFile)

% compiles source data for main figures for sharing in public repo
% each fig has a data structure with panels as fields, and each panel
% contains xaxis, yaxis, and labels as appropriate. Inactivation data
% contain stereotaxic coordinates and pvalues instead, eg:
% these data are deposited in the Mendeley Data platform,
% DOI: 10.17632/d2dkk9647b.1
%
% Fig 1, panel C (performance boxplots)
% Fig1.panelC.visguided   
% Fig1.panelC.memguided   
% Fig1.panelC.accumtowers
% Fig1.panelC.ylabel      
%
% Fig 1, panel G (inactivation maps)
% Fig1.panelG.gridCoord   
% Fig1.panelG.effectSize  
% Fig1.panelG.pvals       

%% FIG 1
load(summaryFile,'taskComp','taskCompStats')
Fig1.panelC.visguided   = taskComp.visGuide.percCorrect.mousePC';
Fig1.panelC.memguided   = taskComp.memGuide.percCorrect.mousePC';
Fig1.panelC.accumtowers = taskComp.towers.percCorrect.mousePC'; 
Fig1.panelC.ylabel      = 'Perf. (% correct)';

Fig1.panelD.visguided   = taskComp.visGuide.viewAngSD.viewAngSD';
Fig1.panelD.memguided   = taskComp.memGuide.viewAngSD.viewAngSD';
Fig1.panelD.accumtowers = taskComp.towers.viewAngSD.viewAngSD'; 
Fig1.panelD.ylabel      = 'View angle SD (deg)';

Fig1.panelF.gridCoord   = taskCompStats.gridCoord;
Fig1.panelF.effectSize  = taskCompStats.visGuideMaze.perf.effSize;
Fig1.panelF.pvals       = taskCompStats.visGuideMaze.perf.pvals;
Fig1.panelF.dataLabel   = 'Norm. \Delta Perf. (%), Vis.-guided task';

Fig1.panelG.gridCoord   = taskCompStats.gridCoord;
Fig1.panelG.effectSize  = taskCompStats.mainMaze.perf.effSize;
Fig1.panelG.pvals       = taskCompStats.mainMaze.perf.pvals;
Fig1.panelG.dataLabel   = 'Norm. \Delta Perf. (%), Accum.-towers task';

Fig1.panelH.gridCoord   = taskCompStats.gridCoord;
Fig1.panelH.effectSize  = taskCompStats.memMaze.perf.effSize;
Fig1.panelH.pvals       = taskCompStats.memMaze.perf.pvals;
Fig1.panelF.dataLabel   = 'Norm. \Delta Perf. (%), Mem.-guided task';

Fig1.panelI_left.visguided   = taskCompStats.normEffectSize.visGuide;
Fig1.panelI_left.memguided   = taskCompStats.normEffectSize.memGuide;
Fig1.panelI_left.accumtowers = taskCompStats.normEffectSize.accumul;
Fig1.panelI_left.ylabel      = 'Norm. laser on perf. (%)';

accumul_pc  = Fig1.panelC.accumtowers;  
visGuide_pc = Fig1.panelC.visguided;  
memGuide_pc = Fig1.panelC.memguided;  
pc          = [nanmean(accumul_pc)*ones(size(Fig1.panelI_left.accumtowers)) nanmean(visGuide_pc)*ones(size(Fig1.panelI_left.visguided)) nanmean(memGuide_pc)*ones(size(Fig1.panelI_left.memguided))]';
lsr         = [Fig1.panelI_left.accumtowers Fig1.panelI_left.visguided Fig1.panelI_left.memguided]';

Fig1.panelI_right.xaxis   = pc;
Fig1.panelI_right.yaxis   = lsr;
Fig1.panelI_right.xlabel  = 'Overall control perf. (% correct)';
Fig1.panelI_right.ylabel  = 'Norm. laser on perf. (%)';

Fig1.note = 'gridCoord is in the form [ML AP] from bregma';

save(figDataFile,'Fig1')
clear taskComp taskCompStats

%% FIG 2
load(summaryFile,'fig2stats')
fieldls = {'panelA_left','panelA_right','panelB','panelC','panelD','panelE','panelF','panelG','panelH','panelI','panelJ'};
for iField = 1:numel(fieldls)
  Fig2.(fieldls{iField}) = fig2stats.(fieldls{iField});
end
save(figDataFile,'Fig2','-append')
clear fig2stats

%% FIG 3
load(summaryFile,'fig3stats')
Fig3 = fig3stats;
save(figDataFile,'Fig3','-append')
clear fig3stats

%% FIG 4
load(summaryFile,'fig4stats')
Fig4 = fig4stats;
save(figDataFile,'Fig4','-append')
clear fig4stats

%% FIG 5
load(summaryFile,'fig5stats')
Fig5 = fig5stats;
save(figDataFile,'Fig5','-append')
clear fig5stats

%% FIG 6
load(summaryFile,'fig6stats')
Fig6 = fig6stats;
save(figDataFile,'Fig6','-append')
clear fig6stats

%% FIG 7
load(summaryFile,'fig7stats')
Fig7 = fig7stats;
save(figDataFile,'Fig7','-append')
clear fig7stats
