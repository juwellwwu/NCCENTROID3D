clc; fclose('all'); clearvars;
close all hidden;

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path
addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Entropy-based leakage image cleanup

% == Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes .MAT file and its corresponding 8-bit image stack)
% "EntropyClean":  2/3D image stack of cleaned vessel leakage


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: June-2021


%% USER INPUT

% === INPUT: Vessel Leakage Image Stack
% InputImg Tif vs MAT Query
InputImg_TifvsMAT_Query=2; % "1" for TIF; "2" for MAT

% Input image directory (do NOT include / at end)
BatchImgInputFolder=''; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
LeakVEOut_Folder_struct=dir(fullfile('OUTPUT/','VesselLeak*'));
ImgStackData_FilenameString=strcat('OUTPUT/',LeakVEOut_Folder_struct(end).name,'/VesselLeak/LeakVE.mat');


% === INPUT: Bone Mask
% Apply bone mask query
% Bone Mask in this script is only used to determine the best z plane for calculating FinalThreshLevel
BoneMask_Query='y';

InputBoneMask_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% Input bone image directory (do NOT include / at end)
% Images insided should be thresholded such that non zero areas = bone
BatchBoneMaskInputFolder='INPUT/DEMO_BoneMask'; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
BoneMaskStackData_FilenameString='';


% === INPUT: Vessel Mask
% Apply vessel mask query
% Vessel Mask in this script is only used to determine the best z plane for calculating FinalThreshLevel 
VesselMask_Query='y';

InputVesselMask_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% Input vessel image directory (do NOT include / at end)
% Images insided should be thresholded such that non zero areas = bone
BatchVesselMaskInputFolder='INPUT/DEMO_VesselMask'; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
VesselMaskStackData_FilenameString='';


% === INPUT: Pixel length in XY, in Z (um)
XY_PxLength=0.3405;
Z_PxLength=1.0;


% === INPUT: Other parameters
% = Output folder string (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/VesselLkClean_Out'; 


% === INPUT: Cell size (um, diameter)
% = Min. cell size (um, diameter)
Cell_diameter_um_min=2.0;  


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

% === INPUT: Entropy Filter specific
% = Z depth in um for Entropy Analysis
Z_Entropy_um=20;

% = Exclusion Query: if 'y', will keep pixels with low entropy but high signal intensity 
Exclusion_Query='n';

% = Exclusion kernel ratio: Kernel size as ratio to Kernel_pix_min 
% Range <= 1.0
% (Entropy Filter size = Kernel_pix_min *1);
Exclusion_Kernel_Ratio=0.5;

% === INPUT: Color of ImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=2;


%% Calculate px size for kernel size

Cell_diameter_px_min=Cell_diameter_um_min/XY_PxLength;
Cell_circumference_px_min=round(2*pi()*Cell_diameter_px_min/2)*0.5; % Full circumference (*1.0)  threshold > half circumference (*0.5)

% Kernel size should be ~1/2 of cell diameter, w/ variance
Kernel_px_min=floor(Cell_diameter_px_min/2-(Cell_diameter_px_min/2/4));


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
else
    ImgMaskStack=ones(size(ImgStack));  % For AllMaskStack; Include all regions for analysis
end

ImgMaskStack(ImgMaskStack<max(ImgMaskStack(:)))=0;
ImgMaskStack=logical(ImgMaskStack);

clearvars BoneMask_Height BoneMask_Width NumBoneMaskSlices;
clearvars VesselMask_Height VesselMask_Width NumVesselMaskSlices;
clearvars ImgMask_Height ImgMask_Width NumImgMaskSlices;


%% Create All Mask Stack: combine info from BoneMaskStack, VesselMaskStack, ImgMaskStack
% AllMaskStack=1 at regions that contain non-bone, non-vessel info 

if (BoneMask_Query=='y') || (VesselMask_Query=='y') || (ImgMask_Query=='y')
    AllMask_Query='y';    
else
    AllMask_Query='n';
end

AllMaskStack=logical(ImgMaskStack) & (~logical(BoneMaskStack)) & (~logical(VesselMaskStack));

clearvars BoneMaskStack VesselMaskStack;


%% Prepare data save directories

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/'); % include "/"
mkdir(SaveFilePath);

fprintf('Preparing save directories...\n');

% SaveEFBWFilePath=strcat(SaveFilePath,'EntropyFilter/BW/');   
% mkdir(SaveEFBWFilePath);
% SaveEFRGBFilePath=strcat(SaveFilePath,'EntropyFilter/RGB/');   
% mkdir(SaveEFRGBFilePath);
% % SaveEFThFilePath=strcat(SaveFilePath,'EntropyThresh/');   
% % mkdir(SaveEFThFilePath);
SaveEFCleanFilePath=strcat(SaveFilePath,'EntropyClean/');   
mkdir(SaveEFCleanFilePath);
% SaveEFCleanRPFilePath=strcat(SaveFilePath,'EntropyCleanRemvPx/');   
% mkdir(SaveEFCleanRPFilePath);


%% Remove Bright Outliers

fprintf('Noise removal of full ImgStack: bright + dark outlier removal...\n');

ImgStack_Input=ImgStack;

ImgStack_NR=zeros(size(ImgStack_Input),'single');

NR_SE=strel('square',3);    % Smallest kernel   
MedianOrd=round(numel(find(NR_SE.Neighborhood))*0.5); % rank of median in neighborhood (~0,5*size in px)

parfor k=1:size(ImgStack_Input,3)
    Img_NR=ImgStack_Input(:,:,k);
    Img_NR_Med=ordfilt2(Img_NR,MedianOrd,NR_SE.Neighborhood);  % median
    Img_NR_Stdev=stdfilt(Img_NR,NR_SE.Neighborhood);   % standard deviation
    Img_NR_Idx=union(find(Img_NR-Img_NR_Med>Img_NR_Stdev*3),find(Img_NR_Med-Img_NR>Img_NR_Stdev*3));  % Hampel Filter (gentle): Px intensity > median + 3*stdev
    % Img_NR_Idx=union(find(Img_NR-Img_NR_Med>0.001),find(Img_NR_Med-Img_NR>0.001));    %ImgJ-style Noise Removal
    Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
    ImgStack_NR(:,:,k)=single(Img_NR);
end

clearvars Img_NR*;
clearvars *_Input;


%% Entropy Filter: Full Stack

fprintf('Running entropy filter on full ImgStack...\n');

ImgStack_Input=ImgStack_NR;

EF_SE=strel('disk',Kernel_px_min);
ImgStack_Entropy=entropyfilt(double(ImgStack_Input),EF_SE.Neighborhood);
ImgStack_Entropy=ImgStack_Entropy./max(ImgStack_Entropy(:));

% Save as BW Image
% figure;
% for k=1:size(ImgStack_Input,3)
%     % imshow(ImgStack_Entropy(:,:,k));
%     imwrite(double(ImgStack_Entropy(:,:,k)),strcat(SaveEFBWFilePath,'EF BW -',num2str(k,'%04.0f'),'.tif'));
% end

% Save as RGB Indexed Image
[HistCt,HistEdge,HistBinIdx] = histcounts(ImgStack_Entropy);
HistBinLoc=0.5*(HistEdge(1:end-1)+HistEdge(2:end));
plot(HistBinLoc(2:end),HistCt(2:end));

cmap=colormap(jet(length(HistEdge)));

ImgStack_Entropy_HistBinRGBIdx=cell(size(ImgStack_Input,3),1);
% for k=1:size(ImgStack_Input,3)
%     ImgStack_Entropy_HistBinRGBIdx{k,1}=ind2rgb(HistBinIdx(:,:,k),cmap);
%     % imshow(EF_HistBinRGBIdx{k,1});
%     imwrite(double(ImgStack_Entropy_HistBinRGBIdx{k,1}),strcat(SaveEFRGBFilePath,'EF RGB -',num2str(k,'%04.0f'),'.tif'));
% end

close all;
pause(1);

clearvars ImgStack_NR;
clearvars ImgStack_Entropy_HistBinRGBIdx;
clearvars Hist* cmap;
clearvars *_Input;


%% Select 21 Z-planes with most pixel intensity information (high occupancy; decent intensity) for entropy analysis

fprintf('Select intensity rich planes from full ImgStack for entropy analysis...\n');

ImgStack_Input=ImgStack;
AllMaskStack_Input=AllMaskStack;
ImgStack_Entropy_Input=ImgStack_Entropy;

ImgStack_Mask0=ImgStack_Input.*AllMaskStack_Input;

NumZplane_Entropy=round(Z_Entropy_um/Z_PxLength);

SignalPxFrac=zeros(NumImgSlices,2);
for k=1:NumImgSlices
    SignalPxFrac(k,1)=k;    % Slice Number
    SignalPxFrac(k,2)=numel(find(ImgStack_Mask0(:,:,k)>0.25))/(size(ImgStack_Input,1)*size(ImgStack_Input,2)-numel(find(~AllMaskStack_Input(:,:,k))));
end

% CRITERIA 1: SignalPxFrac has to be within the values of 10-100% the value
% of its max. This excludes planes of zeros from analysis
[SignalPxFrac_Sort,SignalPxFrac_SortIdx]=sortrows(SignalPxFrac,-2);
SignalPxFrac_Sort_CutoffIdx=min([NumZplane_Entropy,max(find(SignalPxFrac_Sort(:,2)>max(SignalPxFrac(:,2))*0.10)),NumImgSlices]); 
SignalPxFrac_Sort_Trunc=SignalPxFrac_Sort(1:SignalPxFrac_Sort_CutoffIdx,:);

% SignalPxFrac_Trunc_SortIdx holds the ranking of SignalPxFrac for
% each z plane (column 1) in ImgStack. Lower number = higher
% SignalPxFrac (more data rich for entropy analysis)
[SignalPxFrac_Trunc,SignalPxFrac_Trunc_SortIdx]=sortrows(SignalPxFrac_Sort_Trunc,1);
SignalPxFrac_Trunc_k_Frac_FracRank=horzcat(SignalPxFrac_Trunc,SignalPxFrac_Trunc_SortIdx);

% For Saving: Similar to SignalPxFrac_Trunc_k_Frac_FracRank, 
% but include Frac data from all planes
% Column 1: z plane index
% Column 2: SignalPxFrac value
% Column 3: Rank of SignalPxFrac value, smaller = larger Frac (more info)
[~,FracRank]=sortrows(sortrows(SignalPxFrac,-2),1);
SignalPxFrac_k_Frac_FracRank=horzcat(SignalPxFrac,FracRank);

% Truncate
ImgStack_Trunc=ImgStack_Input(:,:,SignalPxFrac_Trunc_k_Frac_FracRank(:,1));
AllMaskStack_Trunc=AllMaskStack_Input(:,:,SignalPxFrac_Trunc_k_Frac_FracRank(:,1));
ImgStack_Entropy_Trunc=ImgStack_Entropy_Input(:,:, SignalPxFrac_Trunc_k_Frac_FracRank(:,1));

clearvars FracRank;
clearvars ImgStack_Mask0;
clearvars SignalPxFrac* -except SignalPxFrac_k_Frac_FracRank SignalPxFrac_Trunc_k_Frac_FracRank;
clearvars SignalPxFrac_Sort_Trunc;
clearvars SignalPxFrac;
clearvars SignalPxFrac_Sort SignalPxFrac_SortIdx SignalPxFrac_Sort_CutoffIdx;
clearvars SignalPxFrac_Trunc SignalPxFrac_Trunc_SortIdx;
clearvars Trunc_MinZ Trunc_MaxZ;
clearvars SignalPxFrac_MovSum;
clearvars *_Input;


%% Threshold of Entropy Filter Image; Loop to find best Threshold

fprintf('Calculating best threshold for entropy filter...\n');

ImgStack_Input=ImgStack_Trunc;
AllMaskStack_Input=AllMaskStack_Trunc;
ImgStack_Entropy_Input=ImgStack_Entropy_Trunc;    % Same size as ImgStack_Input

ImgStack_Entropy_PreThresh=ImgStack_Entropy_Input;
ImgStack_Entropy_PreThresh=ImgStack_Entropy_PreThresh.*AllMaskStack_Input;
ImgStack_Entropy_PreThresh=reshape(ImgStack_Entropy_PreThresh,size(ImgStack_Input,1)*size(ImgStack_Input,2)*size(ImgStack_Input,3),1,1);
ImgStack_Entropy_PreThresh(ImgStack_Entropy_PreThresh==0)=[]; % Ignore all 0 values

ThreshLevel_NoiseFrac=zeros(19,2);

for TLCt=1:39   % Do not use parfor; use parfor at inner level
    
    
        %% Thresholding ImgStack Entropy 
        
        fprintf(strcat('...Thresholding entropy filtered images, iteration ->',num2str(TLCt),'...\n'));
        
        % === Initial Thresholding
        
        % Set Threshold Level
        % [ThreshLevel,~]=graythresh(ImgStack_Entropy_PreThresh);
        ThreshLevel=0.025*TLCt;
        
        ImgStack_EntropyThresh=zeros(size(ImgStack_Entropy_Input),'logical');
        parfor k=1:size(ImgStack_EntropyThresh,3)
            ImgStack_EntropyThresh(:,:,k)=imbinarize(ImgStack_Entropy_Input(:,:,k),ThreshLevel);
        end

        % === Create exclusion pixel list: pixels with low entropy but high brightness are kept
        % High brightness: top 20% intensity of non-zero pixels in ImageStack
        HistCtNumBins=numel(find(ImgStack_Input>0))/1E2; % Ignore 0 pixels
        HistCtEdge=1E-9:max(ImgStack_Input(:))/HistCtNumBins:max(ImgStack_Input(:)); 
        HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
        [HistCt,~]=histcounts(ImgStack,HistCtEdge);

        Hist_SatFrac=20*0.01; % Percentage * 0.01
        [~,SatIndex] = min(abs(cumsum(HistCt)-(1-Hist_SatFrac)*sum(HistCt)));
        SatVal=HistCt_BinLoc(SatIndex);

        % === Final Threshold: Initial Threshold, without Exclusion List
        if Exclusion_Query=='y'
            
            % METHOD 2: Aggressive
            % Consider pixels w/ low entropy and their neighbourhood used for entropy
            % calculations; to avoid edges of leak being included using half size
            % neighbourhood
            Exclusion_SE=strel('disk',round(Kernel_px_min*Exclusion_Kernel_Ratio)); 
            ImgStack_EntropyThresh2=imerode(ImgStack_EntropyThresh,Exclusion_SE.Neighborhood);
            ExclusionIdx=find((ImgStack_Input.*~ImgStack_EntropyThresh2)>SatVal);
            ImgStack_EntropyThreshFinal=ImgStack_EntropyThresh2;
               
            ImgStack_EntropyThreshFinal(ExclusionIdx)=1;            
        
        else
            
            ImgStack_EntropyThreshFinal=ImgStack_EntropyThresh;
            
        end

        clearvars ImgStack_Entropy_PreThresh ImgStack_EntropyThresh2;
        clearvars Hist*;
        clearvars ExclusionIdx;


        %% Clean up based on Thresholded ImgStack Entropy

        fprintf('...Cleaning up, based on thresholded entropy filtered images...\n');

        ImgStack_Input=ImgStack_Trunc;
        ImgStack_EntropyThreshFinal_Input=ImgStack_EntropyThreshFinal;  % Same size as ImgStack_Input

        ImgStack_EntropyClean=single(ImgStack_Input.*ImgStack_EntropyThreshFinal);
        
        
        %% Mapping entropy clean removed pixels. 
        % Removed pixels maps are then analysed for connectivity
        % Final threshold level for full ImgStack is set to avoid removing edges of cells and other features.

        fprintf('...Mapping entropy filter removed pixels...\n');

        ImgStack_Input=ImgStack_Trunc;
        ImgStack_EntropyClean_Input=ImgStack_EntropyClean; % Same size as ImgStack_Input

        % Use bwconncomp to determine size of each "patch" of entropy clean removed
        % pixels. To avoid removing useful edges of cells, the size should ideally be <
        % number of pixels included in the half-circumference of the minimal cell size.
        ImgStack_EntropyCleanRmvPxPreT=abs(ImgStack_Input-ImgStack_EntropyClean_Input);
        ImgStack_EntropyCleanRmvPx_ConnComp=cell(size(ImgStack_Input,3),1);
        parfor k=1:size(ImgStack_EntropyCleanRmvPxPreT,3)
            ImgStack_EntropyCleanRmvPx(:,:,k)=imbinarize(ImgStack_EntropyCleanRmvPxPreT(:,:,k),1E-9);
            ImgStack_EntropyCleanRmvPx_ConnComp{k,1}=bwconncomp(ImgStack_EntropyCleanRmvPx(:,:,k));
        end

        % ConnPxCt_Stack lists the pixel count of each object (=patch of removed pixels).
        % EntropyCleanRmvPx_NoiseFrac reports the fraction of these counts that are
        % smaller than the minimum half-circumference count. This parameter is used to
        % determine the optimal ThreshLevel.
        ConnPxCt_Stack=[];
        parfor k=1:size(ImgStack_Input,3) 
            ConnPxCt_Slice=[];
            for i=1:ImgStack_EntropyCleanRmvPx_ConnComp{k,1}.NumObjects % 1 to NumObjects        
                ConnPxCt_Slice=vertcat(ConnPxCt_Slice,size(ImgStack_EntropyCleanRmvPx_ConnComp{k,1}.PixelIdxList{1,i},1));
            end
            ConnPxCt_Stack=vertcat(ConnPxCt_Stack,ConnPxCt_Slice);
        end

        HistCtEdge=0:1:max(ConnPxCt_Stack(:)); 
        HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
        [HistCt,~]=histcounts(ConnPxCt_Stack,HistCtEdge);

        EntropyCleanRmvPx_NoiseFrac=sum(HistCt(find(HistCt_BinLoc<=Cell_circumference_px_min))',1)/sum(HistCt',1);
        
        clearvars ImgStack_EntropyCleanRmvPxPreT ImgStack_EntropyCleanRmvPx_ConnComp;
        clearvars ConnPxCt_Stack ConnPxCt_Slice;   % Must not be active
        clearvars Hist*;

        ThreshLevel_NoiseFrac(TLCt,:)=[ThreshLevel,EntropyCleanRmvPx_NoiseFrac];

end

clearvars ImgStack_Trunc;
clearvars AllMaskStack_Trunc;
clearvars ImgStack_Entropy_Trunc;
clearvars TLCt;
clearvars *_Input;


%% Determine FInal Threshold Level for full ImgStack

 fprintf('Calculating entropy threshold to be used on full ImgStack...\n');

% Plot: ThreshLevel as X-axis; NoiseFrac as Y-axis
% figure;
% plot(ThreshLevel_NoiseFrac(:,1),ThreshLevel_NoiseFrac(:,2));

[~,NoiseFracMin_RowIdx]=min(ThreshLevel_NoiseFrac(:,2));

if NoiseFracMin_RowIdx==1
    ThreshLevel_NoiseFrac(find(ThreshLevel_NoiseFrac(:,2)<=min(ThreshLevel_NoiseFrac(:,2))+eps),2)=1;
    [~,NoiseFracMin_RowIdx]=min(ThreshLevel_NoiseFrac(:,2));
end

NoiseFracLowerLim=0.95;
NoiseFracLowerLim_RowIdx=find(round(ThreshLevel_NoiseFrac(:,2),2)>=NoiseFracLowerLim);
NoiseFracLowerLim_RowIdx(NoiseFracLowerLim_RowIdx>=NoiseFracMin_RowIdx)=[];

FinalThreshLevel=ThreshLevel_NoiseFrac(max(NoiseFracLowerLim_RowIdx(:)),1);

if isempty(FinalThreshLevel) % Use largest NoiseFrac value
    [~,NoiseFrac_RowIdx]=max(ThreshLevel_NoiseFrac(1:NoiseFracMin_RowIdx,2));
    NoiseFrac_RowIdx=find(ThreshLevel_NoiseFrac(:,2)==ThreshLevel_NoiseFrac(NoiseFrac_RowIdx,2));
    FinalThreshLevel=mean(ThreshLevel_NoiseFrac(NoiseFrac_RowIdx,1),1);
end

clearvars NoiseFracLowerLim*;
clearvars NoiseFrac_RowIdx;


%% Thresholding Entropy-Filtered Full ImgStack

fprintf('Thresholding full ImgStack with pre-calculated entropy threshold...\n');

ImgStack_Input=ImgStack;
ImgStack_Entropy_Input=ImgStack_Entropy;    % Same size as ImgStack_Input

ThreshLevel=FinalThreshLevel;

ImgStack_Entropy_PreThresh=reshape(ImgStack_Entropy_Input,size(ImgStack_Input,1)*size(ImgStack_Input,2)*size(ImgStack_Input,3),1,1);
ImgStack_Entropy_PreThresh(ImgStack_Entropy_PreThresh==0)=[]; % Ignore all 0 values

ImgStack_EntropyThresh=zeros(size(ImgStack_Entropy_Input),'logical');   
for k=1:size(ImgStack_EntropyThresh,3)   % for faster than parfor
    ImgStack_EntropyThresh(:,:,k)=imbinarize(ImgStack_Entropy_Input(:,:,k),ThreshLevel);
end

% === Create exclusion pixel list: pixels with low entropy but high brightness are kept
% High brightness: top 20% intensity of non-zero pixels in ImageStack
HistCtNumBins=numel(find(ImgStack_Input>0))/1E2; % Ignore 0 pixels
HistCtEdge=1E-9:max(ImgStack_Input(:))/HistCtNumBins:max(ImgStack_Input(:)); 
HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
[HistCt,~]=histcounts(ImgStack,HistCtEdge);

Hist_SatFrac=1*0.20; % Percentage * 0.01
[~,SatIndex] = min(abs(cumsum(HistCt)-(1-Hist_SatFrac)*sum(HistCt)));
SatVal=HistCt_BinLoc(SatIndex);

% === Final Threshold: Initial Threshold, without Exclusion List
if Exclusion_Query=='y'
    
    % METHOD 2: Aggressive
    % Consider pixels w/ low entropy and their neighbourhood used for entropy
    % calculations; to avoid edges of leak being included using half size
    % neighbourhood
    Exclusion_SE=strel('disk',round(Kernel_px_min*Exclusion_Kernel_Ratio)); 
    ImgStack_EntropyThresh2=imerode(ImgStack_EntropyThresh,Exclusion_SE.Neighborhood);
    ExclusionIdx=find((ImgStack_Input.*~ImgStack_EntropyThresh2)>SatVal);
    ImgStack_EntropyThreshFinal=ImgStack_EntropyThresh2;

    ImgStack_EntropyThreshFinal(ExclusionIdx)=1;

else
    
    ImgStack_EntropyThreshFinal=ImgStack_EntropyThresh;
    
end

clearvars ImgStack_Entropy;
clearvars ImgStack_Entropy_PreThresh ImgStack_EntropyThresh2;
clearvars ImgStack_EntropyThresh;
clearvars Hist*;
clearvars SatIndex SatVal;
clearvars ExclusionIdx;
clearvars *_Input;


%% Clean up full ImgStack with Entropy Filter

fprintf('Cleanup of full ImgStack, based on thresholded entropy filtered images...\n');

ImgStack_Input=ImgStack;
ImgStack_EntropyThreshFinal_Input=ImgStack_EntropyThreshFinal;  % Same size as ImgStack_Input

ImgStack_EntropyClean=single(ImgStack_Input.*ImgStack_EntropyThreshFinal.*ImgMaskStack);

clearvars ImgMaskStack;
clearvars *_Input;


 %% Mapping entropy clean removed pixels in full ImgStack
 
fprintf('Mapping entropy clean removed pixels in full ImgStack...\n');

ImgStack_Input=ImgStack;
ImgStack_EntropyClean_Input=ImgStack_EntropyClean; % Same size as ImgStack_Input

ImgStack_EntropyCleanRmvPxPreT=abs(ImgStack_Input-ImgStack_EntropyClean_Input);
for k=1:size(ImgStack_EntropyCleanRmvPxPreT,3)  % for faster than parfor
    ImgStack_EntropyCleanRmvPx(:,:,k)=imbinarize(ImgStack_EntropyCleanRmvPxPreT(:,:,k),1E-9);
end

clearvars ImgStack;
clearvars ImgStack_EntropyCleanRmvPxPreT;
clearvars *_Input;


%% Determine % pixels are Removed Pixels, after AllMask 

fprintf('Calculate percent pixels removed by Entropy Filter (after AllMask)...\n');

ImgStack_EntropyCleanRmvPx_Input=ImgStack_EntropyCleanRmvPx;

ImgStack_RmvPxCtFrac=zeros(size(ImgStack_EntropyCleanRmvPx_Input,3),1);
for k=1:size(ImgStack_EntropyCleanRmvPx_Input,3)  % for faster than parfor
    ImgStack_RmvPxCtFrac(k,1)=numel(find(ImgStack_EntropyCleanRmvPx(:,:,k).*logical(AllMaskStack(:,:,k))))/numel(find(AllMaskStack(:,:,k)));
end

clearvars AllMaskStack;
clearvars *_Input;


%% Save: Image Stacks as TIFF and MAT; Workspace

fprintf('Saving results...\n');

% Save final, thresholded Entropy ImgStack
% figure;
% for k=1:size(ImgStack_EntropyThreshFinal,3)
%     % imshow(ImgStack_EntropyThreshFinal(:,:,k));
%     imwrite(double(ImgStack_EntropyThreshFinal(:,:,k)),strcat(SaveEFThFilePath,'EFThreshFinal',num2str(ThreshLevel),' -',num2str(k,'%04.0f'),'.tif'));
% end

% Save cleaned up full ImgStack 
figure;
for k=1:size(ImgStack_EntropyClean,3)
    imshow(ImgStack_EntropyClean(:,:,k).*2);
    imwrite(double(ImgStack_EntropyClean(:,:,k)),strcat(SaveEFCleanFilePath,'EFClean -',num2str(k,'%04.0f'),'.tif'));
end

% Save map of removed pixels in full ImgStack 
% figure;
% for k=1:size(ImgStack_EntropyCleanRmvPx,3)
%     % imshow(ImgStack_EntropyCleanRmvPx(:,:,k));
%     imwrite(double(ImgStack_EntropyCleanRmvPx(:,:,k)),strcat(SaveEFCleanRPFilePath,'EFCleanRmvPx -',num2str(k,'%04.0f'),'.tif'));
% end

close all;
pause(1);

% Save cleaned up full ImgStack as .MAT
ImgStack_EntropyClean=255.*ImgStack_EntropyClean;
save(strcat(SaveEFCleanFilePath,'EntropyClean.mat'),'ImgStack_EntropyClean');

% Save Workspace
save(strcat(SaveFilePath,strcat('Workspace',num2str(timestamp),'.mat')),'SignalPxFrac_k_Frac_FracRank','SignalPxFrac_Trunc_k_Frac_FracRank',...
    'ThreshLevel_NoiseFrac','FinalThreshLevel','ImgStack_RmvPxCtFrac','NumZplane_Entropy','Exclusion_Query','Exclusion_Kernel_Ratio',...
    'AllMask_Query','Cell_diameter_um_min','Cell_circumference_px_min','Kernel_px_min');
% save(strcat(SaveFilePath,'Workspace.mat'),'ThreshLevel_NoiseFrac','SignalPxFrac_MovSum_MaxIdx','FinalThreshLevel','AllMask_Query','Cell_diameter_um_min','Cell_circumference_px_min','Kernel_px_min');

clearvars Batch*Folder;
clearvars Save*Path -except SaveFilePath;
clearvars ImgStack_EntropyThreshFinal;
clearvars ImgStack_EntropyClean;
clearvars ImgStack_EntropyCleanRmvPx;
clearvars *_Input;


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);
