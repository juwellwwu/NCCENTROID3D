clc; fclose('all'); clearvars;
close all hidden;

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path
addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Vessel channel background intensity characterization. 
% Calculates mean and spread, creates mask of regions with insufficient signal for cell counting. 
% Background pixels assumed to be pixels occupied by osteocyte-free bone. 

% == Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes .MAT file and its corresponding 8-bit image stack)
% "BkgdCutoff": 2/3D mask of regions with insufficient signal for cell counting
% "ImgStack_BCAdj": 2/3D image stack of vessels after removing background

% ==Output files:
% "BkgdIntensityEst.mat": background estimate results


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


% === INPUT: Bone Mask
% Apply bone mask query
BoneMask_Query='y';

InputBoneMask_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% Input bone image directory (do NOT include / at end)
% Images insided should be thresholded such that non zero areas = bone
BatchBoneMaskInputFolder='INPUT/DEMO_BoneMask'; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
BoneMaskStackData_FilenameString='';


% === INPUT: Vessel Mask
% Apply vessel mask query
VesselMask_Query='y';

InputVesselMask_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% Input vessel image directory (do NOT include / at end)
% Images insided should be thresholded such that non zero areas = bone
BatchVesselMaskInputFolder='INPUT/DEMO_VesselMask'; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
VesselMaskStackData_FilenameString='';


% === INPUT: Output folder
% Output folder string (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/BkgdEstCutoff_Out'; 


% === Pixel length in XY, in Z (um)
XY_PxLength=0.3405; 
Z_PxLength=1.0;


% === INPUT: Cell size (um, diameter)
% = Min. cell size (um, diameter)
Cell_diameter_um_min=2.0;  

% = Max. cell size (um, diameter)
Cell_diameter_um_max=20;



%% Optimized inputs: modify with care

% === INPUT: Img Mask (Rotation)
ImgMask_Query='n';

% Input image mask directory (do NOT include / at end)
InputImgMask_TifvsMAT_Query=1;

% Input image mask directory (do NOT include / at end)
% Image mask = 1 where signal is to be kept
BatchImgMaskInputFolder=''; 

% Input mask images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
ImgMaskStackData_FilenameString='';

% === INPUT: Img Mask (Region)
% Single Slice, restrict XY region for cell count
% 1 for sub-region to count; 0 to ignore
RegionMask_Query='n';

% Input image mask directory (do NOT include / at end)
InputRegionMask_TifvsMAT_Query=1;

% Input image mask directory (do NOT include / at end)
% Image mask = 1 where signal is to be kept
BatchRegionMaskInputFolder=''; 

% Input mask images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
RegionMaskStackData_FilenameString='';

% === INPUT: Manually Background Intensity Estimate
% Used to calculate ThreshCutoff if no bone mask supplied, or if no pixel available for 
% bone-mask-determined ThreshCutOff
% If maximum intensity in local neighbourhood (determined by cell diameter) 
% < ThreshCutoff=Background+ThreshCutoff_Padding, set pixel to 0
% NOTE: ThreshCutoff=ManualIntensityCutoff+ThreshCutoff_Padding
ManualBkgdEst=0.062; % Get from ImageJ
ManualBkgdStdevEst=0.0051; % Get from ImageJ

% === INPUT: Color of ImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=2;


%% Determine and verify stack file count

% Load image stack
[ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad_f(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel);

% Load bone mask stack, if needed
if BoneMask_Query=='y'
    [BoneMaskStack,BoneMask_Height,BoneMask_Width,NumBoneMaskSlices]=ImgStackLoad_f(InputBoneMask_TifvsMAT_Query,BatchBoneMaskInputFolder,BoneMaskStackData_FilenameString,RGBChannel);
    if (NumBoneMaskSlices~=NumImgSlices) || (BoneMask_Height~=Img_Height) || (BoneMask_Width~=Img_Width)
        error('Error. BoneMaskStack size must equal ImgStack size.')
    end
else
    BoneMaskStack=zeros(size(ImgStack));    % For AllMaskStack; Assume no bone
end

% Load vessel mask stack, if needed
if VesselMask_Query=='y'
    [VesselMaskStack,VesselMask_Height,VesselMask_Width,NumVesselMaskSlices]=ImgStackLoad_f(InputVesselMask_TifvsMAT_Query,BatchVesselMaskInputFolder,VesselMaskStackData_FilenameString,RGBChannel);
    if (NumVesselMaskSlices~=NumImgSlices) || (VesselMask_Height~=Img_Height) || (VesselMask_Width~=Img_Width)
        error('Error. VesselMaskStack size must equal ImgStack size.')
    end
else
    VesselMaskStack=zeros(size(ImgStack));   % For AllMaskStack; Assume no vessel
end

% Load image mask stack, if needed
if ImgMask_Query=='y'
    [ImgMaskStack,ImgMask_Height,ImgMask_Width,NumImgMaskSlices]=ImgStackLoad_f(InputImgMask_TifvsMAT_Query,BatchImgMaskInputFolder,ImgMaskStackData_FilenameString,RGBChannel);
    if (NumImgMaskSlices~=NumImgSlices) || (ImgMask_Height~=Img_Height) || (ImgMask_Width~=Img_Width)
        error('Error. ImgMaskStack size must equal ImgStack size.')
    end
    ImgMaskStack=logical(ImgMaskStack);
else
    ImgMaskStack=ones(size(ImgStack),'logical');  % For AllMaskStack; Include all regions for analysis
end

% Load region mask (single slice), if needed
if RegionMask_Query=='y'
    [RegionMaskStack,RegionMask_Height,RegionMask_Width,NumRegionMaskSlices]=ImgStackLoad_f(InputRegionMask_TifvsMAT_Query,BatchRegionMaskInputFolder,RegionMaskStackData_FilenameString,RGBChannel);
    if (NumRegionMaskSlices~=1) || (RegionMask_Height~=Img_Height) || (RegionMask_Width~=Img_Width)
        error('Error. RegionMask must be single Z slice; [Rol,Col] must equal ImgStack [Row,Col].')
    end
    RegionMaskStack=logical(repmat(RegionMaskStack,1,1,NumImgSlices));
else
    RegionMaskStack=logical(repmat(ones(Img_Height,Img_Width),1,1,NumImgSlices));  % For AllMaskStack; Include all regions for analysis
end

clearvars *_Height *_Width Num*Slices -except Img_Height Img_Width NumImgSlices;


%% Create save directories. Use timestamp in name to prevent overwrite

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/'); % include "/"
mkdir(SaveFilePath);


%% Remove Bright /Dim Outliers

fprintf('Noise removal of full ImgStack: bright + dark outlier removal...\n');

ImgStack_Input=ImgStack;

ImgStack_NR=zeros(size(ImgStack_Input),'single');

NR_SE=strel('disk',round(Cell_diameter_um_min*0.5/XY_PxLength));    % Smallest kernel   
MedianOrd=round(numel(find(NR_SE.Neighborhood))*0.5); % rank of median in neighborhood (~0,5*size in px)

parfor k=1:NumImgSlices
    Img_NR=ImgStack_Input(:,:,k);
    Img_NR_Med=ordfilt2(Img_NR,MedianOrd,NR_SE.Neighborhood);  % median
    Img_NR_Stdev=stdfilt(Img_NR,NR_SE.Neighborhood);   % standard deviation
    Img_NR_Idx=union(find(Img_NR-Img_NR_Med>Img_NR_Stdev*3),find(Img_NR_Med-Img_NR>Img_NR_Stdev*3));  % Hampel Filter (gentle): Px intensity > median + 3*stdev
    % Img_NR_Idx=union(find(Img_NR-Img_NR_Med>0.001),find(Img_NR_Med-Img_NR>0.001));    %ImgJ-style Noise Removal
    Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
    ImgStack_NR(:,:,k)=single(Img_NR);
end

clearvars Img_NR*;
clearvars MedianOrd;
clearvars *_Input;


%% Estimate Background Value, using Bone Mask
% Get background estimate from bone regions

ImgStack_Input=ImgStack_NR;
BoneMaskStack_Input=BoneMaskStack;

if BoneMask_Query=='y' 
    
    % Erode Bone Mask by ~ multiples of cell width to avoid segmentation imperfections
    % Use 2D: rough est sufficient, imerode works with 3D but very slow
    Erode_SE=strel('disk',round(5*0.5/XY_PxLength)*3);
    BoneMaskStack_Erode=zeros(size(BoneMaskStack_Input),'logical');
    for k=1:NumImgSlices
        BoneMaskStack_Erode(:,:,k)=imerode(BoneMaskStack_Input(:,:,k),Erode_SE);
    end
    
    % Bone Criteria: Consider leakage channel signal in bone regions only
    BoneMaskStack_ZMaxStack=logical(repmat(min(BoneMaskStack_Erode,[],3),1,1,NumImgSlices)); % Stack of ZMax of BoneMask
    
    % ImgStack Criteria: avoid osteocytes
    ImgStack_ZMax=max(ImgStack_Input./max(ImgStack_Input(:)),[],3); % Stack of ZMax of ImgStack
    ImgStack_ZMax(ImgStack_ZMax>0.10)=NaN; % Avoid osteocyte regions; assume background no larger than 0.10
    ImgStack_ZMaxStack=repmat(ImgStack_ZMax,1,1,NumImgSlices);
    
    % Estimate background = mean intensity in leakage channel in bone
    % regions w/out osteocytes
    ImgStack_BkgdEst=ImgStack_Input./max(ImgStack_Input(:));
    ImgStack_BkgdEst(~ImgMaskStack)=NaN;
    ImgStack_BkgdEst(~BoneMaskStack_ZMaxStack)=NaN;
    ImgStack_BkgdEst(find(isnan(ImgStack_ZMaxStack)))=NaN;
    
    figure;
    for k=1:NumImgSlices
        imshow(double(ImgStack_BkgdEst(:,:,k)*10)); 
    end

    ImgStack_BkgdEst=reshape(ImgStack_BkgdEst,numel(ImgStack_BkgdEst),1,1);
    ImgStack_BkgdEst(find(isnan(ImgStack_BkgdEst)))=[];
    ImgStack_BkgdEst(ImgStack_BkgdEst<eps)=[];
    
    if numel(ImgStack_BkgdEst)>100
    
        % Histogram: distribution of intensities
        [HistCt,HistEdge]=histcounts(ImgStack_BkgdEst,200);
        HistBinLoc=0.5*(HistEdge(1:end-1)+HistEdge(2:end));

        HistData=[HistCt',HistBinLoc'];
        HistData=HistData(find(HistData(:,1)),:); % Eliminate any zero rows from histogram stretching
        
        % Try to fit; if not, use manual values
        try
            
            HistFitObject=fit(double(HistData(:,2)),double(HistData(:,1)),'gauss1');
            
%             figHandle01=figure;
%             plot(HistFitObject,HistData(:,2),HistData(:,1));
%             ylabel('Histogram Ct'); 
%             xlabel('Pixel Intensity');
%             title({strcat('Vessel Channel Pixel Intensity in Bone Region')});
%             print(figHandle01,'-dtiffn','-r0',strcat(SaveFilePath,'VesselChannelIntensityHistogram (BoneRegionOnly).tif'));

            % Gaussian centroid, sigma, FWHM
            BkgdCentroid=HistFitObject.b1;
            BkgdSigma=sqrt(HistFitObject.c1^2/2);
            BkgdFWHM=2*sqrt(2*log(2))*BkgdSigma;
            
        catch
            
            fprintf('...Gaussian Fit cannot be completed, likely due to low background. Use Mean and Stdev estimate... \n');
            HistFitObject=[];
            ImgStack_BkgdEst(ImgStack_BkgdEst>0.25)=[];
            BkgdCentroid=mean(ImgStack_BkgdEst(:));
            BkgdSigma=std(ImgStack_BkgdEst(:));
            BkgdFWHM=2*sqrt(2*log(2))*BkgdSigma;
            
        end
    
        fprintf('\n=====================\n');
        fprintf(strcat('Centroid of Bkgd Gaussian Fit=',num2str(BkgdCentroid),' (in 8-bit scale=',num2str(BkgdCentroid*255),')\n'));
        fprintf(strcat('Sigma of Bkgd Gaussian Fit=',num2str(BkgdSigma),' (in 8-bit scale=',num2str(BkgdSigma*255),')\n'));
        fprintf(strcat('FWHM of Bkgd Gaussian Fit=',num2str(BkgdFWHM),' (in 8-bit scale=',num2str(BkgdFWHM*255),')\n'));
        fprintf('=====================\n\n');
        
    else
        
        HistFitObject=[];
        BkgdCentroid=ManualBkgdEst;
        BkgdSigma=ManualBkgdStdevEst;
        BkgdFWHM=2*sqrt(2*log(2))*BkgdSigma;
        
    end
    
else
    
    HistFitObject=[];
    BkgdCentroid=ManualBkgdEst;
    BkgdSigma=ManualBkgdStdevEst;
    BkgdFWHM=2*sqrt(2*log(2))*BkgdSigma;
        
end

clearvars BoneMaskStack_Erode;
clearvars BoneMaskStack_ZMaxStack;
clearvars ImgStack_ZMax ImgStack_ZMaxStack;
clearvars ImgStack_BkgdEst;
clearvars *Hist* -except HistFitObject;
clearvars *_Input;


%% Calculate Background Intensity 

fprintf('Estimating Background intensity; at least 10 in 8-bit scale...\n');

BkgdEst=(BkgdCentroid+BkgdFWHM/2)+(10/255); % (10/255) as "padding"

fprintf('\n=====================\n');
fprintf(strcat('Bkgd Intensity Estimate=',num2str(BkgdEst),' (in 8-bit scale=',num2str(BkgdEst*255),')\n'));
fprintf('=====================\n\n');


%% Create list of background pixels
 
ImgStack_Input=ImgStack;

if ImgMask_Query=='y'
    ImgStack_Input(~ImgMaskStack)=9999;
end

ImgStack_BkgdLinIdxList=find(ImgStack_Input<BkgdEst);

clearvars *_Input;


%% Save background estimate information

save(strcat(SaveFilePath,'BkgdIntensityEst.mat'),...
    'HistFitObject','BkgdEst','BkgdCentroid','BkgdSigma','BkgdFWHM','ImgStack_BkgdLinIdxList');

clearvars HistFitObject;
clearvars ImgStack_BkgdLinIdxList;


%% Creating Intensity Cutoff Mask

fprintf('Creating Intensity Cutoff Mask...\n');

ImgStack_Input=ImgStack_NR;

% Apply masks
if ImgMask_Query=='y'
    ImgStack_Input=ImgStack_Input.*ImgMaskStack;
end

if RegionMask_Query=='y'
    ImgStack_Input=ImgStack_Input.*RegionMaskStack;
end

% Minimally Gaussian blur image
ImgStack_Gauss=zeros(size(ImgStack));
for k=1:NumImgSlices
     ImgStack_Gauss(:,:,k)=imgaussfilt(ImgStack_Input(:,:,k),Cell_diameter_um_min/XY_PxLength/6);
end

% Dilate vascular mask as vessel signal high shifts maximum intensity
HalfMaxDiameter_SE=strel('disk',round(Cell_diameter_um_max*0.5*0.5/XY_PxLength));

if VesselMask_Query=='y'
    VesselMaskStack_MaxFilt=logical((imdilate(VesselMaskStack,HalfMaxDiameter_SE.Neighborhood)));
else
    VesselMaskStack_MaxFilt=zeros(Img_Height,Img_Width,NumImgSlices,'logical');
end

ImgStack_PreMaxFilt=ImgStack_Gauss;
ImgStack_PreMaxFilt(VesselMaskStack_MaxFilt)=0;
ImgStack_MaxFilt=single(imdilate(ImgStack_PreMaxFilt,HalfMaxDiameter_SE.Neighborhood));

% Eliminate Low signal
ThreshCutoff=BkgdEst;
ImgStack_MaxFiltThresh=zeros(size(ImgStack_MaxFilt),'logical');
ImgStack_MaxFiltThresh(ImgStack_MaxFilt>ThreshCutoff-eps)=1;

% Create Background Intensity Cutoff Mask
ImgStack_BCMask=logical(ImgStack_MaxFiltThresh+VesselMaskStack_MaxFilt);

% Apply masks again
if ImgMask_Query=='y'
    ImgStack_BCMask=ImgStack_BCMask.*ImgMaskStack;
end

if RegionMask_Query=='y'
    ImgStack_BCMask=ImgStack_BCMask.*RegionMaskStack;
end

clearvars ImgStack_NR ImgStack_Gauss;
clearvars VesselMaskStack_MaxFilt;
clearvars ImgStack_MaxFilt* ImgStack_PreMaxFilt;
clearvars RegionMaskStack VesselMaskStack;
clearvars *_Input;


%% Save Intensity Cutoff Mask

fprintf('Saving Background Cutoff Mask...\n');

SaveBCMaskFilePath=strcat(SaveFilePath,'BkgdCutoff',num2str(ThreshCutoff),'/');
mkdir(SaveBCMaskFilePath);

for k=1:NumImgSlices
    imwrite(double(ImgStack_BCMask(:,:,k)*255),strcat(SaveBCMaskFilePath,'BCMask -',num2str(k,'%04.0f'),'.tif'));
end

Save_ImgStack_BCMask=ImgStack_BCMask.*255;
save(strcat(SaveBCMaskFilePath,'BCMask.mat'),'Save_ImgStack_BCMask');

clearvars ImgStack_BCMask Save_ImgStack_BCMask;


%% Brightness Contrast Adjustment, Original Image Stack

ImgStack_Input=ImgStack;

ImgStack_Input=ImgStack_Input./max(ImgStack_Input(:));
ImgStack_BCAdj=zeros(size(ImgStack_Input),'single');
for k=1:NumImgSlices
    ImgStack_BCAdj(:,:,k)=imadjust(ImgStack_Input(:,:,k),[BkgdEst 1],[0 1]);
end

clearvars ImgStack_Input;


%% Save Brightness Contrast Adjusted Original Image Stack

fprintf('Saving Brightness Contrast Adjusted Original Image Stack...\n');

SaveBCAdjFilePath=strcat(SaveFilePath,'ImgStack_BCAdj/');
mkdir(SaveBCAdjFilePath);

for k=1:NumImgSlices
    imwrite(double(ImgStack_BCAdj(:,:,k)),strcat(SaveBCAdjFilePath,'BCAdj -',num2str(k,'%04.0f'),'.tif'));
end

Save_ImgStack_BCAdj=ImgStack_BCAdj.*255;
save(strcat(SaveBCAdjFilePath,'BCAdj.mat'),'Save_ImgStack_BCAdj');

clearvars ImgStack_BCAdj Save_ImgStack_BCAdj;
clearvars figHandle*;
clearvars k timestamp;

clearvars ImgStack ImgMaskStack BoneMaskStack VesselMaskStack;


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);

