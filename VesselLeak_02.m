clc; fclose('all'); clearvars;
close all hidden; 

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path
addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Leakage image creation. Subtract VesselEst from vessel channel images. 

% == Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes .MAT file and its corresponding 8-bit image stack)
% "VesselLeak":  2/3D image stack of vessel leakage 


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: June-2021


%% USER INPUT


% === INPUT: Original Vessel Channel Image Stack
% InputImg Tif vs MAT Query:
InputImg_TifvsMAT_Query=2; % "1" for TIF; "2" for MAT

% Input image directory (do NOT include / at end)
BatchImgInputFolder=''; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
ImgStackData_FilenameString='INPUT/DEMO_ImgStack_G.mat';


% === INPUT: VesselEst
% InputImg Tif vs MAT Query:
InputVesselEst_TifvsMAT_Query=2; % "1" for TIF; "2" for MAT

% Input image directory (do NOT include / at end)
BatchImgVesselEstInputFolder=''; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
VesselEstOut_Folder_struct=dir(fullfile('OUTPUT/','VesselEst*'));
ImgVesselEstStackData_FilenameString=strcat('OUTPUT/',VesselEstOut_Folder_struct(end).name,'/VesselEst/VesselEst.mat');

% === Output image directory (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/VesselLeak_Out'; 


% === INPUT: Pixel length in XY, in Z (um)
XY_PxLength=0.3405;
Z_PxLength=1.0;


% === INPUT: Average feature size (um, diameter)
Diameter_um_av=4.0; % Larger than largest leak, smaller than sinusoids and transitionals


%% Optimized inputs: modify with care

% === INPUT: Img Mask
% Apply image mask query
% Image mask should be "1" where signal is to be kept
ImgMask_Query='n';

% Input image mask directory (do NOT include / at end)
InputImgMask_TifvsMAT_Query=1;

% Input image mask directory (do NOT include / at end)
% Image mask should be white (RGB) where signal is to be kept
BatchImgMaskInputFolder=''; 

% Input mask images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
ImgMaskStackData_FilenameString='';

% === INPUT: Color of ImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=1;


%% Load Input Img Stack and ImgMask Stack using functions

% Load image stack
[ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad_f(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel);

% Load VesselEst stack
[ImgVesselEstStack,ImgVesselEst_Height,ImgVesselEst_Width,NumImgVesselEstSlices]=ImgStackLoad_f(InputVesselEst_TifvsMAT_Query,BatchImgVesselEstInputFolder,ImgVesselEstStackData_FilenameString,RGBChannel);
if (NumImgVesselEstSlices~=NumImgSlices) || (ImgVesselEst_Height~=Img_Height) || (ImgVesselEst_Width~=Img_Width)
        error('Error. ImgVesselEstStack size must equal ImgStack size.')
end

% Load image mask stack, if needed
if ImgMask_Query=='y'
    [ImgMaskStack,ImgMask_Height,ImgMask_Width,NumImgMaskSlices]=ImgStackLoad_f(InputImgMask_TifvsMAT_Query,BatchImgMaskInputFolder,ImgMaskStackData_FilenameString,RGBChannel);
    if (NumImgMaskSlices~=NumImgSlices) || (ImgMask_Height~=Img_Height) || (ImgMask_Width~=Img_Width)
        error('Error. ImgMaskStack size must equal ImgStack size.')
    end
else
    ImgMaskStack=[];
end

ImgMaskStack(ImgMaskStack<max(ImgMaskStack(:)))=0;
ImgMaskStack=logical(ImgMaskStack);

clearvars ImgMask_Height ImgMask_Width NumImgMaskSlices;


%% METHOD 1: Subtract Vessel Estimate Stack from Input ImgStack
% This is similar to an alternate rolling ball background correction

fprintf('VESSELEST: Subtracting VesselEst ...\n');

ImgStack_Input=ImgStack;

ImgStack_VE=single(ImgStack_Input-ImgVesselEstStack);

% Rescale matrix such that lowest val=0; highest val=1
ImgStack_VE_Bkgd=(0-min(ImgStack_VE(:)))*1/(max(ImgStack_VE(:))-min(ImgStack_VE(:)));   % for re-scale at final step
ImgStack_VE=(ImgStack_VE-min(ImgStack_VE(:)))*1/(max(ImgStack_VE(:))-min(ImgStack_VE(:)));

clearvars *_Input;


%% METHOD 1: Mild noise removal (after VesselEst subtraction) to clear intensity spikes

fprintf('VESSELEST: Post Noise Removal (Hampel Filter)...\n');

ImgStack_Input=ImgStack_VE;

% strel_Radius is smaller value allowed: gentle noise removal
strel_Radius=3;  
SE=strel('disk',strel_Radius);
MedianOrd=round(numel(find(SE.Neighborhood))*0.5); % rank of median in neighborhood (~0,5*size in px)

% === Noise Removal
% Use ordfilt2 (instead of medfilt2) to allow neighborhood shape variety;
% Neighborhood should be symmetric in this case
ImgStack_NR=zeros(size(ImgStack_Input));
parfor k=1:size(ImgStack_Input,3)
    Img_NR=ImgStack_Input(:,:,k);
    Img_NR_Med=ordfilt2(Img_NR,MedianOrd,SE.Neighborhood);  % median
    Img_NR_Stdev=stdfilt(Img_NR,SE.Neighborhood);   % standard deviation
    Img_NR_Idx=union(find(Img_NR-Img_NR_Med>Img_NR_Stdev*3),find(Img_NR_Med-Img_NR>Img_NR_Stdev*3));  % Hampel Filter (gentle): Px intensity > median + 3*stdev
    % Img_NR_Idx=union(find(Img_NR-Img_NR_Med>0.2),find(Img_NR_Med-Img_NR>0.2));   % ImageJ-style Noise Removal: Px intensity > 50+median in uint8 scale  
    Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
    ImgStack_NR(:,:,k)=Img_NR;
end

clearvars ImgStack_VE;
clearvars Img_NR*
clearvars *_Input;


%% METHOD 1: Leakage by VesselEst Subtraction -- Final

ImgStack_Input=ImgStack_NR;

% Set ImgStack_VE_Bkgd as 0 intensity
ImgStack_ZeroSet=single(ImgStack_Input-ImgStack_VE_Bkgd);
ImgStack_ZeroSet(ImgStack_ZeroSet<=0)=0;
ImgStack_ZeroSet=ImgStack_ZeroSet.*1./max(ImgStack_ZeroSet(:)); % Rescale

if ImgMask_Query=='y'
    ImgStack_ZeroSet(~ImgMaskStack)=0;
end

% Histogram determine saturated value for rescale
ImgStack_Input=ImgStack_ZeroSet;

ImgStack_Input_Reshape=reshape(ImgStack_Input,numel(ImgStack_Input),1,1);
ImgStack_Input_Reshape_LowHighIntensity= stretchlim(ImgStack_Input_Reshape,[0 1-0.0001]);
ImgStack_Input_Sat=ImgStack_Input_Reshape_LowHighIntensity(2,1);

ImgStack_SatRescale=single(ImgStack_Input.*1./ImgStack_Input_Sat);
ImgStack_SatRescale(ImgStack_SatRescale>1)=1;

if ImgMask_Query=='y'
    ImgStack_SatRescale(~ImgMaskStack)=0;
end

% Output
ImgStack_Leak_VE=ImgStack_SatRescale;

clearvars ImgStack_NR;
clearvars ImgStack_VE_Bkgd;
clearvars ImgStack_SatRescale;
clearvars ImgStack_ZeroSet;
clearvars ImgStack_Input_Reshape*
clearvars HistCtNumBins HistCtEdge HistCt_BinLoc;
clearvars ImgStack_Input_HistCt;
clearvars ImgStack_Input_Sat;
clearvars *_Input;


%% Display and Save

fprintf('Saving...\n');

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/'); 
mkdir(SaveFilePath);

% === VE
% Create subfolders for rolling ball and rolling ball background
SaveLeakVEFilePath=strcat(SaveFilePath,'VesselLeak/');
mkdir(SaveLeakVEFilePath);

% Save as images
for k=1:NumImgSlices
    imwrite(double(ImgStack_Leak_VE(:,:,k)),strcat(SaveLeakVEFilePath,'LeakVE -',num2str(k,'%04.0f'),'.tif'));
end

% Save as matrices
% Multiply by 255 (value similar to uint8) for easy understanding of numericals
Save_ImgStack_Leak_VE=255.*ImgStack_Leak_VE;
save(strcat(SaveLeakVEFilePath,'LeakVE.mat'),'Save_ImgStack_Leak_VE');

clearvars Save_ImgStack_Leak_VE;
clearvars ImgStack ImgVesselEstStack ImgMaskStack; % Initial user inputs


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);

