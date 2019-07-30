function Pinto2019_plotAllPaperFigs

%%
% runs all analyses and plots all figures going in the optogenetic
% inactivation and widefield paper (Pinto et al, 2019, Task-dependent
% changes and the large-scale dynamics and necessity of cortical regions)
% Code written by Lucas Pinto (lpinto@princeton.edu)
%
% Most dependencies are github and bitbucket repositories owned by either 
% BRAIN CoGS (https://github.com/BrainCOGS), Lucas Pinto (lpinto@princeton.edu)
% or Sue Ann Koay (koay@princeton.edu). The bitbucket ones are private, 
% access must be requested via email for now.
% The dependencies apply mostly for analyzing data from scratch, the main
% dependency for this repository is publicationscripts, which manages
% plotting and version control thereof.
%
%   THIS REPO: https://github.com/BrainCOGS/Pinto2019_plotAllPaperFigs.git
%   BEHAVIORAL AND INACTIVATION ANALYSIS: https://github.com/BrainCOGS/behavioralAnalysis.git
%   IMAGING ANALYSIS: https://github.com/BrainCOGS/widefieldImaging.git
%   LASER RIG CONTROL AND ANALYSIS DEPENDENCIES: https://github.com/BrainCOGS/laserGalvoControl.git
%   OTHERS:
%   https://sakoay@bitbucket.org/sakoay/publicationscripts.git (figure plotting)
%   https://sakoay@bitbucket.org/sakoay/sakfunctions.git (various)
%   https://sakoay@bitbucket.org/sakoay/tankmousevr.git (VR and behavior dependencies)
%   https://sakoay@bitbucket.org/sakoay/tankmouseanalysis.git (various)
%   https://github.com/BaselLaserMouse/AllenBrainAPI.git (Allen API interface)
% 
% Note that in the behavior and imaging repos the user should edit the class
% objects analysisParams and widefieldParams, respectively,
% such that paths reflect the user's machine. 
%
% Other dependencies for data collection / generation / preprocessing,
% but not high level analysis:
%   https://github.com/kwikteam/klusta.git (spike sorting)
%   https://sakoay@bitbucket.org/sakoay/isivstim.git (retinotopic mapping)
%   https://github.com/sakoay/princeton-ecs.git (imaging preprocessing)
%
% Additionally, some matlab file exchange packages are also required:
%   https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
%   https://www.mathworks.com/matlabcentral/fileexchange/23342-real2rgb---colormaps
%   https://www.mathworks.com/matlabcentral/fileexchange/27812-rotatexlabels--ax--angle--varargin--/content/rotateXLabels.m
%
% The functions used for analysis are listed inside each figure-specific function,
% but realistically they must be run on the cluster (defaults are for PNI's spock), 
% running everything locally would take weeks. All relevant functions are 
% compatible with spock and most are parallelized. 
% Workflow for low-level widefield data preprocessing is provided in
% Pinto2019_fig3_WFdynamics.m

%% Set paths
% If you are at PNI:
% summary files are already in braininit bucket, but you also have the
% option of loading from local disk for speed
% change paths for your local directory structure
loadFromBucket     = true; % set to true if loading from bucket server
if loadFromBucket
  analysisPathOpto = '/Volumes/braininit/Analysis/laserGalvo/';
  analysisPathWF   = '/Volumes/braininit/RigData/VRwidefield/widefield/';
  analysisPathRNN  = '/Volumes/braininit/RigData/Analysis/RNN/';
  ephysPath        = '/Volumes/braininit/RigData/ephys/';
  varPowerPath     = '/Volumes/braininit/Analysis/laserGalvo/V1varPower.mat';
else
  analysisPathOpto = []; 
  analysisPathWF   = [];
  analysisPathRNN  = [];
  ephysPath        = [];
  varPowerPath     = [];
end

% location of .mat file that stores summary analysis results. If empty, the
% summary file will not be saved
summaryFile        = '/Volumes/braininit/Analysis/Pinto2019/Pinto2019_analysis.mat';
figDataFile        = '/Volumes/braininit/Analysis/Pinto2019/Pinto2019_figSourceData.mat';

%% MAIN FIG: schematics + behav summary + whole-trial inactivation 
Pinto2019_fig1_wholeTrialOpto(analysisPathOpto,summaryFile);

%% supplemental: comparison of behavioral performance between the tasks
Pinto2019_figS1_taskPerfComp(analysisPathOpto,summaryFile);
Pinto2019_figS1_speedSingleTrialEg(analysisPathOpto);

%% supplemental: opto characterization ephys, var power, WT controls
Pinto2019_figS2_ephys_varPower(ephysPath,varPowerPath,analysisPathOpto,summaryFile)
Pinto2019_figS2_wholeTrialOpto_WTctrl(analysisPathOpto);

%% MAIN FIG: other indicators and clustering, towers and mem guided task
Pinto2019_fig2_wholeTrialOptoClust(analysisPathOpto,summaryFile);

%% supplemental: subtrial cue period + more inactivation maps (whole trial), combined location psychometrics etc
Pinto2019_figS3_subTrialOpto(analysisPathOpto,summaryFile);
Pinto2019_figS3_wholeTrialOpto_extraInfo(analysisPathOpto,summaryFile);

%% MAIN FIG: wf dynamics  (accum-towers vs vis-guided: time progression) 
Pinto2019_fig3_WFdynamics(analysisPathWF,summaryFile);

%% supplemental: hemodynamic correction (method, airpuff, yfp)
Pinto2019_figS4_hemodynamicCorrection(analysisPathWF,summaryFile);

%% supplemental: ROI selection, comp. w/ maps 
Pinto2019_figS5_ROI(analysisPathWF,summaryFile);

%% supplemental: accum-towers vs vis-guided: time progression scatter plots, pixel "sequences"
Pinto2019_figS6_WFdynamics_epochComp(analysisPathWF); % 1st half
Pinto2019_figS6_WFpxlSeq(analysisPathWF, summaryFile); % 2nd half
 
%% MAIN FIG: wf correlations  (accum-towers vs vis-guided: corr) 
Pinto2019_fig4_WFcorr(analysisPathWF,summaryFile);

%% supplemental: wf correlations vs. distance, finer scales etc
Pinto2019_figS7_WFcorr(analysisPathWF,summaryFile);

%% MAIN FIG: evidence tuning
Pinto2019_fig5_evidenceTuning(analysisPathWF,summaryFile);

%% MAIN FIG: decoding
Pinto2019_fig6_WFdecoding(analysisPathWF,summaryFile);

%% supplemental: view angle regression control + task comp
Pinto2019_figS8_WFdecodingViewAngCtrl(analysisPathWF,summaryFile);

%% supplemental: triggered averages + GLM
Pinto2019_figS9_WFtriggeredAvg(analysisPathWF,summaryFile);
Pinto2019_figS9_WFglm(analysisPathWF,summaryFile);

%% MAIN FIG: RNN
Pinto2019_fig7_RNN(analysisPathRNN,summaryFile);

%% supplemental: RNN
Pinto2019_figS10_RNN_extraInfo(analysisPathRNN,summaryFile);

%% generate figure source data (main figs only)
Pinto2019_compileFigSourceData(figDataFile,summaryFile)

