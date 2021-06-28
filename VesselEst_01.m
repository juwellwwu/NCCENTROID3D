clc; fclose('all'); clearvars;
close all hidden;

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Vessel estimate: to be subtracted from the vessel channel images to get the leakage images.

% == Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes .MAT file and its corresponding 8-bit image stack)
% "VesselEst": 2/3D image stack of vessel estimate 


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: June-2021


%% USER INPUT

% === INPUT: Original Vessel Channel Image Stack
% InputImg Tif vs MAT Query:
InputImg_TifvsMAT_Query=2;

% Input image directory (do NOT include / at end)
BatchImgInputFolder='';

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255)  
ImgStackData_FilenameString='INPUT/DEMO_ImgStack_G.mat'; 


% === INPUT: Output image directory (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/VesselEst_Out'; 


% === INPUT: Pixel length in XY, in Z (um)
XY_PxLength=0.3405; 
Z_PxLength=1.0;

% Min. and average feature size (um, diameter)
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
RGBChannel=2;


%% Load Input Img Stack and ImgMask Stack using functions

% Load image stack
[ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad_f(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel);

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


%% Create Crude (Conservative) Vessel Estimation by Median-Min-Max Filter Seq
% Step1: NoiseReduction

fprintf('Estimating Vessels w/ NoiseReduction-Median-Min-Max Sequence...\n');
fprintf('Step1: NoiseReduction...\n');

ImgStack_Input=ImgStack;

% Create neighborhood (must be flat) for all operations
% strel_Radius should reflect estimated largest leakage < smaller vessel size
strel_size=round(Diameter_um_av/XY_PxLength);
if mod(strel_size,2)==0 % strel size must be odd
    strel_size=strel_size-1;
end
SE=strel('square',strel_size);
MedianOrd=round(numel(find(SE.Neighborhood))*0.5); % rank of median in neighborhood (~0,5*size in px)

% === Noise Removal
% Use ordfilt2 (instead of medfilt2) to allow neighborhood shape variety;
% Neighborhood should be symmetric in this case
ImgStack_NR=zeros(size(ImgStack_Input));

tic
parfor k=1:NumImgSlices
    Img_NR=ImgStack_Input(:,:,k);
    Img_NR_Med=ordfilt2(Img_NR,MedianOrd,SE.Neighborhood);  % median
    Img_NR_Stdev=stdfilt(Img_NR,SE.Neighborhood);   % standard deviation
    % Img_NR_Idx=find(Img_NR-Img_NR_Med>Img_NR_Stdev*3);  % Hampel Filter (gentle): Px intensity > median + 3*stdev
    Img_NR_Idx=find(Img_NR-Img_NR_Med>0.2);   % ImageJ-style Noise Removal: Px intensity > 50+median in uint8 scale  
    Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
    ImgStack_NR(:,:,k)=Img_NR;
end
toc

clearvars ImgStack;
clearvars Img_NR*;
clearvars *_Input;


%% Create Crude (Conservative) Vessel Estimation by Median-Min-Max Filter Seq
% Step2: NoiseReduction

fprintf('Step2: Median...\n');

ImgStack_Input=ImgStack_NR;

% === Median Filtering
% Use ordfilt2 (instead of medfilt2) to allow neighborhood shape variety;
% Neighborhood should be symmetric in this case
ImgStack_Median=zeros(size(ImgStack_Input));

tic
parfor k=1:NumImgSlices   % parfor faster than for
    ImgStack_Median(:,:,k)=ordfilt2(ImgStack_Input(:,:,k),MedianOrd,SE.Neighborhood); 
end
toc

% Reset ImgMask areas to 0
if ImgMask_Query=='y'
    ImgStack_Median(~ImgMaskStack)=0;
end

clearvars ImgStack_NR;
clearvars *_Input;


%% Create Crude (Conservative) Vessel Estimation by Median-Min-Max Filter Seq
% Step3: Minimum filter

fprintf('Step3: Minimum filter (imerode)...\n');

ImgStack_Input=ImgStack_Median;

% === Minimum Filtering
% MATLAB Greyscale minimum: imerode
ImgStack_Input=ImgStack_Median;
ImgStack_Min=zeros(size(ImgStack_Input));

tic
for k=1:NumImgSlices   % for faster than parfor
    ImgStack_Min(:,:,k)=imerode(ImgStack_Input(:,:,k),SE);
end
toc

% Reset ImgMask areas to 0
if ImgMask_Query=='y'
    ImgStack_Min(~ImgMaskStack)=0;
end

clearvars ImgStack_Median;
clearvars *_Input;


%% Create Crude (Conservative) Vessel Estimation by Median-Min-Max Filter Seq
% Step3: Minimum filter

fprintf('Step4: Maximum filter (imdilate)...\n');

ImgStack_Input=ImgStack_Min;

% === Maximum Filtering
% MATLAB Greyscale maximum: imdilate
ImgStack_Max=zeros(size(ImgStack_Input));

tic
for k=1:NumImgSlices   % for faster than parfor
    ImgStack_Max(:,:,k)=imdilate(ImgStack_Input(:,:,k),SE);
end
toc

% Reset ImgMask areas to 0
if ImgMask_Query=='y'
    ImgStack_Max(~ImgMaskStack)=0;
end

% === Final Output
ImgStack_VesselEst=single(ImgStack_Max);

clearvars ImgStack_Min;
clearvars ImgStack_Max;
clearvars ImgMaskStack;
clearvars *_Input;


%% Display and Save
% "VesselEst -XXXX": Vessel Channel, cavity region rough vessel estimate

fprintf('Saving...\n');

timestamp=datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/'); 
mkdir(SaveFilePath);

SaveVesselEstFilePath=strcat(SaveFilePath,'VesselEst/');
mkdir(SaveVesselEstFilePath);

% Save as images
for k=1:size(ImgStack_VesselEst,3)
    imwrite(double(ImgStack_VesselEst(:,:,k)),strcat(SaveVesselEstFilePath,'VesselEst -',num2str(k,'%04.0f'),'.tif'));
end

% Save as matrices
% Multiply by 255 (value similar to uint8) for easy understanding of numericals
Save_ImgStack_VesselEst=255.*ImgStack_VesselEst;

save(strcat(SaveVesselEstFilePath,'VesselEst.mat'),'Save_ImgStack_VesselEst');

clearvars Save_ImgStack_VesselEst;


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);

