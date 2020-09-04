%% Script to calculate the Nuclear to Cytosolic intensity (NCI) of a given signal in images
%%, here labeled YAP*. The images NUC* allow to detect and segment the
%%cell's nuclei

%%The output results are given in a .mat file. 

%%For information about how to cite this material, please go to the
%%Readme.txt file. 




clear; clc;
close all;


%%%OUTPUTFOLDER

pathOutput =pwd;

%%%MATLAB ROUTINES FOLDER

pathCalculations=pwd; 


%PARAMETERS ALLOWING TO SEGMENT THE NUCLEI AND CALCULATE NCI 

sigma = 4;          % blurring for FINE nuclei identification and tracking [px], 
minSize = 1000;   % minimum size for nuclei identification [px^2]
maxSize=15000; %This is also useful because sometimes two cells are attached (SZ)
nclusters = 1;  % Clusters of the image used for the rought segmentation. The threshold is the mean of the centroids of the intensities from the lower to the ncluster-1 max of the value
treshFactor_FINE = 0.5;  % threshold factor for nuclei identification
TolArea=0.25; 

%%%These are for the ring of NCI

RingWidth=30;
factorTolBG=5; 
sigmaforring=1;

parCellTrack = [sigma, minSize, nclusters, treshFactor_FINE, maxSize,TolArea,RingWidth,factorTolBG,sigmaforring]; % Parameters for cell tracking



pathInputQuant=strcat(pwd,'\');
pathInputTrack=strcat(pwd,'\');
cd(pathInputQuant)
fnamesQUANT = dir('YAP*.tif');
cd(pathInputTrack)
fnamesTRACK = dir('NUC*.tif');

numfids = length(fnamesQUANT);
cd(pathCalculations);


for n=1:numfids


fileNameQuant =fnamesQUANT(n).name; 
fileNameTrack=fnamesTRACK(n).name; 



cd(pathCalculations);



filenameprint=strrep(fileNameQuant,'.tif','');
filenamesave=strrep(fileNameQuant,'.tif','.mat');



[StackQuant, nFrames] = TIFread([pathInputQuant, fileNameQuant]);

[StackTrack, nFrames] = TIFread([pathInputTrack, fileNameTrack]);

whattoprint='BOTHANDRINGS'; 

Track_Quant_for_NCI_Ring_classic_WS_broad(StackTrack,StackQuant, parCellTrack,filenameprint, filenamesave,whattoprint,pathOutput);



end; 
