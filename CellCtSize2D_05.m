clc; fclose('all'); clearvars;
close all hidden; 

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path
addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Watershed-based 2D cell counting. 
% Provides information on cell density, cell diameter, bone, vessel and mask areas for each Z slice.

% == Usage: 
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes an 8-bit image stack)
% "Colocal WSLines OrgImgEC": Co-localized watershed lines and original vessel channel image stack
% "Colocal WSMap OrgImgEC": Co-localized watershed map and original vessel
% channel image stack. Watershed map includes watershed lines and coloring by bone / vessel / masks.
% "Colocal WSCentroid OrgImgEC": Co-localized watershed centroids and original vessel channel image stack
% "GdMap": Watershed map re-drawn in voxel-size consistent grid; input for CellCtSize3D. 

% ==Output files:
% "CellCtAreaDensityPxUnit AllSliceData.txt": 2D cell count, area, density
% data in pixel unit
% "CellCtAreaDensityMetric AllSliceData.txt": 2D cell count, area, density
% data in metric unit
% "WatershedCellPosArea AllSliceData.txt": 2D cell area and centroid position data
% "WSCellSizeHistogram.tif": 2D cell size distribution histogram


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: June-2021


%% USER INPUT

% === INPUT: Vessel Leakage Image Stack 
% Input Img Tif vs MAT Query:
InputImg_TifvsMAT_Query=2; % "1" for TIF; "2" for MAT

% Input image directory (do NOT include / at end)
BatchImgInputFolder=''; 

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
VesselLkCleanOut_Folder_struct=dir(fullfile('OUTPUT/','VesselLkClean*'));
ImgStackData_FilenameString=strcat('OUTPUT/',VesselLkCleanOut_Folder_struct(end).name,'/EntropyClean/EntropyClean.mat');

% === INPUT: Original Vessel Channel Image Stack
% Input OrgImg Tif vs MAT Query:
InputOrgImg_TifvsMAT_Query=2;

% Input original image directory (do NOT include / at end)
BatchOrgImgInputFolder=''; 

% Input original images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
OrgImgStackData_FilenameString='INPUT/DEMO_ImgStack_G.mat';

% === INPUT: Output folder
% Output folder string (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/CellCtSize2D_Out';

% === INPUT: Bone Mask
% Apply bone mask query
BoneMask_Query='y';

InputBoneMask_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% Input bone image directory (do NOT include / at end)
% Images insided should be thresholded such that non zero areas = bone
BatchBoneMaskInputFolder='INPUT/DEMO_BoneEndostealMask'; 

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


% === INPUT: Brightness-Contrast Mask (Low Intensity Removal)
% Bcakground Intensity Cutoff Mask: Only regions with sufficiently high intensity for
% analysis specified
% 1 for sub-region to count; 0 to ignore
BCMask_Query='y';

% Input image mask directory (do NOT include / at end)
InputBCMask_TifvsMAT_Query=2;

% Input image mask directory (do NOT include / at end)
% Image mask = 1 where signal is to be kept
BatchBCMaskInputFolder=''; 

% Input mask images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255) 
BkgdEstCutoffOut_Folder_struct=dir(fullfile('OUTPUT/','BkgdEstCutoff*'));
BCMaskStackOut_Folder_struct=dir(fullfile(strcat('OUTPUT/',BkgdEstCutoffOut_Folder_struct(end).name),'BkgdCutoff*'));
BCMaskStackData_FilenameString=...
    strcat('OUTPUT/',BkgdEstCutoffOut_Folder_struct(end).name,'/',BCMaskStackOut_Folder_struct(end).name,'/BCMask.mat');


% === INPUT: Background Intensity Value
% Used for cleaning up VesselLkClean stack
% Values used for creating BCMask
BC_Clean_Query='y';
BC_BkgdIntensityEst_FilenameString=strcat('OUTPUT/',BkgdEstCutoffOut_Folder_struct(end).name,'/BkgdIntensityEst.mat'); 

% === INPUT: Pixel length in XY, in Z (um)
XY_PxLength=0.3405;
Z_PxLength=1.0;


% === INPUT: CellCt3D_Query
% If 'y', data is prepared for 3D Cell Counting
% * Set Cell_diameter_um_min=size equivalent to smallest kernel size;
% * Relax watershed property flagging @ regionprops
CellCt3D_Query='y';

% === INPUT: Cell size (um, diameter)
% = Min. cell size (um, diameter)
% ** For 2D cell ct, do not go lower than 2um **
% Check for Gaussian shape in WS Cell Size Histogram
% Value ignored if CellCt3D_Query=='y'
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

% % === INPUT: Use ImageJ's "Find Maxima" to find centroids
% (Not optimized)
FindMaxima_Query='n';

% === INPUT: Color of ImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=2;

% === INPUT: Structural Element
% = Flat vs Non-Flat Query
% Non-Flat gives higher cts, but not necessary accurate
FlatStrel_Query='y';    
% strel_Radius=3;     % Radius of flat element (disk)
% offsetstrel_Radius=3;       % Radius of non-flat element (ball), >3 for effect 
offsetstrel_MaxOffset=0.001;  % Max offset of non-flat element (ball), larger offset->image brighter


%% Calculate px size for kernel size, counting

if CellCt3D_Query=='y'
    Cell_diameter_px_min=3; % minimal kernel size limit
else
    Cell_diameter_px_min=Cell_diameter_um_min/XY_PxLength;
end
Cell_diameter_px_max=Cell_diameter_um_max/XY_PxLength;

% Kernel size should be ~1/2 of cell diameter, w/ variance 
% OR 3, whichever is greater
Kernel_px_min=max(3,floor(Cell_diameter_px_min/2-(Cell_diameter_px_min/2/4)));
Kernel_px_max=max(3,ceil(Cell_diameter_px_max/2+(Cell_diameter_px_min/2/4)));

% Px area for analyze particles
Px_Area_min=round((Cell_diameter_px_min/2).^2*pi);
Px_Area_max=round((Cell_diameter_px_max/2).^2*pi);

fprintf(strcat('Min kernel size suggested for cell ct (particle analysis)-> ',num2str(Kernel_px_min),'\n'));
fprintf(strcat('Max kernel size suggested for cell ct (particle analysis)-> ',num2str(Kernel_px_max),'\n'));
fprintf(strcat('Note: Min kernel size is limited to >=3px, corresponding to cell diameter (um)-> ',num2str(4+(Cell_diameter_px_min/2/4)*2),'\n\n'));

fprintf(strcat('Min pixel area suggested for cell ct (particle analysis)-> ',num2str(Px_Area_min),'\n'));
fprintf(strcat('Max pixel area suggested for cell ct (particle analysis)-> ',num2str(Px_Area_max),'\n\n'));


%% Determine and verify stack file count

% Load image stack
[ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad_f(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel);

% Load orginal image stack 
[OrgImgStack,OrgImg_Height,OrgImg_Width,NumOrgImgSlices]=ImgStackLoad_f(InputOrgImg_TifvsMAT_Query,BatchOrgImgInputFolder,OrgImgStackData_FilenameString,RGBChannel);
if (NumOrgImgSlices~=NumImgSlices) || (OrgImg_Height~=Img_Height) || (OrgImg_Width~=Img_Width)
    error('Error. OrgImgStack size must equal ImgStack size.')
end

% Load bone mask stack, if needed
if BoneMask_Query=='y'
    [BoneMaskStack,BoneMask_Height,BoneMask_Width,NumBoneMaskSlices]=ImgStackLoad_f(InputBoneMask_TifvsMAT_Query,BatchBoneMaskInputFolder,BoneMaskStackData_FilenameString,RGBChannel);
    if (NumBoneMaskSlices~=NumImgSlices) || (BoneMask_Height~=Img_Height) || (BoneMask_Width~=Img_Width)
        error('Error. BoneMaskStack size must equal ImgStack size.')
    end
    BoneMaskStack=logical(BoneMaskStack);
else
    BoneMaskStack=zeros(size(ImgStack),'logical');    % For AllMaskStack; Assume no bone
end

% Load vessel mask stack, if needed
if VesselMask_Query=='y'
    [VesselMaskStack,VesselMask_Height,VesselMask_Width,NumVesselMaskSlices]=ImgStackLoad_f(InputVesselMask_TifvsMAT_Query,BatchVesselMaskInputFolder,VesselMaskStackData_FilenameString,RGBChannel);
    if (NumVesselMaskSlices~=NumImgSlices) || (VesselMask_Height~=Img_Height) || (VesselMask_Width~=Img_Width)
        error('Error. VesselMaskStack size must equal ImgStack size.')
    end
    VesselMaskStack=logical(VesselMaskStack);
else
    VesselMaskStack=zeros(size(ImgStack),'logical');   % For AllMaskStack; Assume no vessel
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

% Load intensity cutoff mask stack, if needed
if BCMask_Query=='y'
    [BCMaskStack,BCMask_Height,BCMask_Width,NumBCMaskSlices]=ImgStackLoad_f(InputBCMask_TifvsMAT_Query,BatchBCMaskInputFolder,BCMaskStackData_FilenameString,RGBChannel);
    if (NumBCMaskSlices~=NumImgSlices) || (BCMask_Height~=Img_Height) || (BCMask_Width~=Img_Width)
        error('Error. BCMaskStack size must equal ImgStack size.')
    end
    BCMaskStack=logical(BCMaskStack);
else
    BCMaskStack=ones(size(ImgStack),'logical');  % For AllMaskStack; Include all regions for analysis
end

clearvars *_Height *_Width Num*Slices -except Img_Height Img_Width NumImgSlices;


%% Expand ImgMaskStack to suppress edge artifacts

Erode_SE=strel('disk',Kernel_px_min);
ImgMaskStack=logical(imerode(ImgMaskStack,Erode_SE));

clearvars Erode_SE;


%% Load Background Intensity Value Estimate Information; 
% Also define BC_BkgdNoiseRange

if BC_Clean_Query=='y'
    
    fprintf('Loading background intensity estimate Information...\n');
    
    BC_BkgdIntensityEst=load(BC_BkgdIntensityEst_FilenameString);
    BC_BkgdEst=BC_BkgdIntensityEst.BkgdEst;
    BC_BkgdCentroid=BC_BkgdIntensityEst.BkgdCentroid;
    BC_BkgdSigma=BC_BkgdIntensityEst.BkgdSigma;
    BC_BkgdFWHM=BC_BkgdIntensityEst.BkgdFWHM;
    BC_BkgdLinIdxList=BC_BkgdIntensityEst.ImgStack_BkgdLinIdxList;

else
    
    BC_BkgdEst=0;
    BC_BkgdCentroid=0;
    BC_BkgdSigma=0;
    BC_BkgdFWHM=0;
    BC_BkgdLinIdxList=[];

end

% BC_BkgdNoiseRange=BC_BkgdFWHM;
BC_BkgdNoiseRange=BC_BkgdSigma*5; % 20181217: BC_BkgdSigma*4 to BC_BkgdSigma*5 optimal


%% Screensize, for image saving. Do Not Remove.

scrsz=get(groot,'ScreenSize');
iptsetpref('ImshowBorder','tight');


%% Prepare results matrix, output directories

% = Create save directories. Use timestamp in name to prevent overwrite

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/'); % include "/"
mkdir(SaveFilePath);

SavePreWSFilePath=strcat(SaveFilePath,'PreWS/');    % Last ImgStack modification before watershed
SaveWSLineFilePath=strcat(SaveFilePath,'Watershed Line/');  % Watershed lines only
SaveWSFgdMFilePath=strcat(SaveFilePath,'WS_FGdMarker/');    % Foreground markers, marker-controlled watershed
SaveWSBgdMFilePath=strcat(SaveFilePath,'WS_BkGdMarker/');  % Background markers, marker-controlled watershed
SaveColocalCellPosImgFilePath=strcat(SaveFilePath,'Colocal WSCentroid Img/');   % FindMaxima, Watershed Centroid Coordinate Colocalization
SaveColocalCellPosOrgImgFilePath=strcat(SaveFilePath,'Colocal WSCentroid OrgImgEC/');
SaveWSMapFilePath=strcat(SaveFilePath,'WSMap/');
SaveColocalWSMapImgFilePath=strcat(SaveFilePath,'Colocal WSMap Img/'); % Cell/Bone/Vessel Mapped Watershed Regions Colocalization
SaveColocalWSMapOrgImgFilePath=strcat(SaveFilePath,'Colocal WSMap OrgImgEC/'); 
SaveGdMapFilePath=strcat(SaveFilePath,'GdMap/');
SaveColocalGdMapImgFilePath=strcat(SaveFilePath,'Colocal GdMap Img/'); % Cell/Bone/Vessel Mapped Grid Regions Colocalization
SaveColocalGdMapOrgImgFilePath=strcat(SaveFilePath,'Colocal GdMap OrgImgEC/'); 
SaveColocalCellBoundImgFilePath=strcat(SaveFilePath,'Colocal WSLine Img/'); % Cell Boundaries Colocalization
SaveColocalCellBoundOrgImgFilePath=strcat(SaveFilePath,'Colocal WSLine OrgImgEC/'); 

% mkdir(SavePreWSFilePath);
% mkdir(SaveWSLineFilePath);
% mkdir(SaveWSFgdMFilePath);
% mkdir(SaveWSBgdMFilePath);
% mkdir(SaveColocalCellPosImgFilePath);
mkdir(SaveColocalCellPosOrgImgFilePath);
% mkdir(SaveWSMapFilePath);
% mkdir(SaveColocalWSMapImgFilePath);
mkdir(SaveColocalWSMapOrgImgFilePath);
mkdir(SaveGdMapFilePath);
% mkdir(SaveColocalGdMapImgFilePath);
% mkdir(SaveColocalGdMapOrgImgFilePath);
% mkdir(SaveColocalCellBoundImgFilePath);
mkdir(SaveColocalCellBoundOrgImgFilePath);


%% Mild Clean up with Gaussian Filtering

fprintf('Cleanup ImgStack: mild gaussian filtering...\n');

ImgStack_Input=ImgStack;

ImgStack_Gauss=zeros(size(ImgStack_Input),'single');

Gauss_sigma=(Cell_diameter_px_min)/6;

for k=1:NumImgSlices
    ImgStack_Gauss(:,:,k)=imgaussfilt(ImgStack_Input(:,:,k),Gauss_sigma);
end

ImgStack_Gauss=ImgStack_Gauss./max(ImgStack_Gauss(:));

clearvars Gauss_sigma;
clearvars _Input;


%% Cleanup by Intensity Range

fprintf('Cleanup ImgStack: by intensity range...\n');

ImgStack_Input=ImgStack_Gauss;

MaxMin_SE=strel('disk',Kernel_px_min);
MedianOrd=round(numel(find(MaxMin_SE.Neighborhood))*0.5); % rank of median in neighborhood (~0,5*size in px)

Img_MaxFilt=zeros(Img_Height,Img_Width,'single');
Img_MinFilt=zeros(Img_Height,Img_Width,'single');
Img_RangeFilt=zeros(Img_Height,Img_Width,'single');

ImgStack_RangeFlagCt=0;
ImgStack_RangeFlagCt_Old=-1;

while (ImgStack_RangeFlagCt-ImgStack_RangeFlagCt_Old)>eps
    
    fprintf('.');
    
    ImgStack_RangeFlagCt_Old=ImgStack_RangeFlagCt;
    ImgStack_RangeFlagCt=0;

    for k=1:NumImgSlices
        
        Img_Input=ImgStack_Input(:,:,k);
        
        Img_MaxFilt=imdilate(Img_Input,MaxMin_SE);
        Img_MinFilt=imerode(Img_Input,MaxMin_SE);
        Img_MedFilt=ordfilt2(Img_Input,MedianOrd,MaxMin_SE.Neighborhood);  % Median

        Img_RangeFilt=Img_MaxFilt-Img_MinFilt;
        Img_RangeFlagID=intersect(find(eps<Img_RangeFilt),find(Img_RangeFilt<BC_BkgdNoiseRange));

        Img=zeros(Img_Height,Img_Width,'logical');
        Img(Img_RangeFlagID)=1;
        Img=logical(imdilate(Img,MaxMin_SE));

        Img_RangeFlagID=find(Img); % Update background pixel list
        
        if BC_Clean_Query=='y'  % Include bkgd pixels in median value replace list
            Img_BC_BkgdLinIdxList=BC_BkgdLinIdxList-Img_Height*Img_Width*k;
            Img_BC_BkgdLinIdxList(Img_BC_BkgdLinIdxList<1)=[];
            Img_BC_BkgdLinIdxList(Img_BC_BkgdLinIdxList>Img_Height*Img_Width)=[];
            Img_RangeFlagID=unique(union(Img_RangeFlagID,Img_BC_BkgdLinIdxList));
        end

        ImgStack_RangeFlagCt=ImgStack_RangeFlagCt+numel(Img_RangeFlagID);
        Img_Input(Img_RangeFlagID)=Img_MedFilt(Img_RangeFlagID);
        ImgStack_Input(:,:,k)=Img_Input;

    end

end

fprintf('Done.\n');
ImgStack_IRClean=ImgStack_Input;

clearvars ImgStack_Gauss;
clearvars Img*Filt;
clearvars ImgStack_RangeFlagID ImgStack_RangeFlagCt*;
clearvars Img Img_Input;
clearvars *_Input;


%% Stack Prep: XY Mean Filter

ImgStack_Input=ImgStack_IRClean;

[ImgStack_XYAv]=ImgStackXYAv_f(ImgStack_Input,'n',[],FlatStrel_Query,Kernel_px_min,Kernel_px_min);

% Rescale matrix such that lowest val=0; highest val=1
ImgStack_XYAv=(ImgStack_XYAv-min(ImgStack_XYAv(:)))*1/(max(ImgStack_XYAv(:))-min(ImgStack_XYAv(:)));

clearvars ImgStack_BCClean;
clearvars ImgStack_IRClean;
clearvars *_Input;


%% Stack Prep: XY Min Filter

ImgStack_Input=ImgStack_XYAv;

[ImgStack_XYMin]=ImgStackXYMin_f(ImgStack_Input,'n',[],FlatStrel_Query,Kernel_px_min,Kernel_px_min,offsetstrel_MaxOffset);

% Rescale matrix such that lowest val=0; highest val=1
ImgStack_XYMin=(ImgStack_XYMin-min(ImgStack_XYMin(:)))*1/(max(ImgStack_XYMin(:))-min(ImgStack_XYMin(:)));

clearvars ImgStack_XYAv;
clearvars *_Input;


%% Combine ImgMaskStack and RegionMaskStack

ImgMaskStack=single(min(ImgMaskStack,RegionMaskStack));

clearvars RegionMaskStack;


%% Stack Prep: Invert LUT (Output: PreWS Stack)

fprintf('Inverting...\n');

ImgStack_Input=ImgStack_XYMin;  % Input should span range 0-1

ImgStack_Invert=single(abs(1-ImgStack_Input));  % Can also use imcomplement()

clearvars ImgStack_XYMin;
clearvars *_Input;


%% Prepare cells and matrices to record watershed results 

% 1) ImgStack of Watershed region labels
% 2) Properties of each Watershed Region in each Z 

% 3) Watershed Regions in each Z disqualified from bone mask only
% 4) Watershed Regions in each Z disqualified from vessel mask only
% 5) Watershed Regions in each Z disqualified from brightness cutoff mask only
% 6) Watershed Regions in each Z disqualified from img mask only

% 7) Watershed Regions in each Z disqualified from regionprops only
% 8) Watershed Regions in each Z disqualified from regionprops + bone/vessel/BC/ImgMask
% 9) Watershed Regions in each Z QUALIFIED as cells

% 10) Centroid of Watershed regions identified as cells, in XY coordinates, after regionprops and
% bone/vessen/ImgMask disqualification
% 11) Centroid of Watershed regions identified as cells, in linear indices, after regionprops and
% bone/vessen/ImgMask disqualification
% 12) Area of Watershed regions identified as cells, in pixel ct, after regionprops and
% bone/vessen/ImgMask disqualification

% 13) Number of Watershed Regions in each Z disqualified from regionprops
% only, after each marker-controlled watershed iteration

ImgStack_WS_Lb=zeros(Img_Height,Img_Width,NumImgSlices);
WS_regionprops_Cell=cell(NumImgSlices,1);

WS_BoneDQ_FlagID_Cell=cell(NumImgSlices,1);
WS_VesselDQ_FlagID_Cell=cell(NumImgSlices,1);
WS_BCMaskDQ_FlagID_Cell=cell(NumImgSlices,1);
WS_ImgMaskDQ_FlagID_Cell=cell(NumImgSlices,1);

WS_PropDQ_FlagID_Cell=cell(NumImgSlices,1);
WS_DQ_FlagID_Cell=cell(NumImgSlices,1);

WS_Cell_FlagID_Cell=cell(NumImgSlices,1);

WS_Cell_Centroid_XY_List_Cell=cell(NumImgSlices,1);
WS_Cell_Centroid_LinIdx_List_Cell=cell(NumImgSlices,1);
WS_Cell_Area_List_Cell=cell(NumImgSlices,1);

WS_Prop_FlagID_Ct_Mtx=nan(NumImgSlices,50);


%% Marker-Controlled Seperation using Watershed 
% Preparation of watershed follow instructions on:
% https://www.mathworks.com/help/images/marker-controlled-watershed-segmentation.html
% Option 1+3:open+close w/out reconstruction doesn't work well; remove

fprintf('Marker-controlled watershed...\n');

ImgStack_Input=ImgStack_Invert;  % Input should span range 0-1

% == Step2: Gradient
ImgStack_PreGradient=zeros(Img_Height,Img_Width,NumImgSlices,'single');
for k=1:NumImgSlices
   ImgStack_PreGradient(:,:,k)= imgradient(ImgStack_Input(:,:,k));
end

if BC_Clean_Query=='y'    
    ImgStack_PreGradient(BC_BkgdLinIdxList)=0;
end

% == Step3: Mark Foreground     
% Option2+4: opening-by-reconstruction (OBR) followed by closing-by-reconstruction (OBR)
% OBR = imerode followed by imreconstruct
% CBR = imdilte followed by imreconstruct

FGd_SE=strel('disk',Kernel_px_min);           
ImgStack_OBRCBR=zeros(Img_Height,Img_Width,NumImgSlices,'single');
for k=1:NumImgSlices
   Img_Erode=imerode(ImgStack_Input(:,:,k),FGd_SE);
   Img_OBR=imreconstruct(Img_Erode,ImgStack_Input(:,:,k));
   Img_OBR_Dilate=imdilate(Img_OBR,FGd_SE);
   ImgStack_OBRCBR(:,:,k)=imcomplement(imreconstruct(imcomplement(Img_OBR_Dilate),imcomplement(Img_OBR)));
end

% == Foreground marker
% Create, Cleanup, Collect List of Marker size
ImgStack_FgdMarker=zeros(Img_Height,Img_Width,NumImgSlices,'single');
FGdMarkerClean_SE=strel('square',3);
FgdMarker_Area_List=[];
for k=1:NumImgSlices                
    Img_FgdMarker=imregionalmax(ImgStack_OBRCBR(:,:,k));    % Create
    Img_Close=imclose(Img_FgdMarker,FGdMarkerClean_SE); % Clean up
    ImgStack_FgdMarker(:,:,k)=logical(imerode(Img_Close,FGdMarkerClean_SE));            
    [Img_FgdMarker_Lb,~]=bwlabel(ImgStack_FgdMarker(:,:,k));    % Collect foreground marker size
    FgdMarker_regionprops=regionprops(Img_FgdMarker_Lb,'Area');
    FgdMarker_Area_List=vertcat(FgdMarker_Area_List,cat(1,FgdMarker_regionprops.Area));           
end            
FgdMarker_Area_List_Min=quantile(FgdMarker_Area_List,0.005);  

% == Step4: Mark Background
ImgStack_BkgdMarker=zeros(Img_Height,Img_Width,NumImgSlices,'single');

for k=1:NumImgSlices
    Img_OBRCBR_Thresh=imbinarize(ImgStack_OBRCBR(:,:,k));   % Threshold
    Img_OBRCBR_DistMap_WS=watershed(bwdist(Img_OBRCBR_Thresh));   % Distance Map & Watershed
    ImgStack_BkgdMarker(:,:,k)=Img_OBRCBR_DistMap_WS==0;    % Background marker
end

% == Step5: Watershed
% Restrict regional minima to only foreground and background markers

% = Redo Gradient (outside iteration)
ImgStack_Gradient=zeros(Img_Height,Img_Width,NumImgSlices,'single');
for k=1:NumImgSlices
   % ImgStack_Gradient(:,:,k)= imimposemin(ImgStack_PreGradient(:,:,k),ImgStack_BkgdMarker(:,:,k)|ImgStack_FgdMarker2(:,:,k));
   ImgStack_Gradient(:,:,k)= imimposemin(ImgStack_PreGradient(:,:,k),ImgStack_FgdMarker(:,:,k));
end

ItCt=1;
WS_It_ZFlagID=(1:1:NumImgSlices)';
WS_It_ZFlagID_EqualCheck=0;
while ItCt>0
        
        % == Step5: Watershed
        % = Watershed (inside iteration)
        ImgStack_WS=zeros(Img_Height,Img_Width,NumImgSlices,'logical');
        for k=1:NumImgSlices
           ImgStack_WS(:,:,k)= watershed(ImgStack_Gradient(:,:,k));
        end

        % = Clean up watershed (not from MATLAB)
        % "imgradient" creates edge-like structure in WS image that
        % affects distance measurement. Cleanup by image open and perform
        % watershed again
        WSClean_SE=strel('square',Kernel_px_min*2);
        for k=1:NumImgSlices
           Img_WS_Erode=imcomplement(imerode(ImgStack_WS(:,:,k),WSClean_SE));
           ImgStack_WS(:,:,k)=watershed(Img_WS_Erode);  % Redo watershed
        end

        % = Overlay Bone, Vessel, Masks Boundaries on Watershed
        % Ensure watershed regions are present for bone and vessels
        ImgStack_WS_BVMEdge=ImgStack_WS;        
        for k=1:NumImgSlices 
            MaskStackEdge=logical(imerode(abs(logical(...
                edge(BoneMaskStack(:,:,k))+...
                edge(VesselMaskStack(:,:,k))+...
                edge(ImgMaskStack(:,:,k)))-1),strel('square',3)));
            ImgStack_WS_BVMEdge(:,:,k)=min(ImgStack_WS_BVMEdge(:,:,k),MaskStackEdge);
        end

        
      %% Cell ID from Regions identified by Watershed in MATLAB 

        fprintf(strcat('\nAnalyzing watershed cells: Iteration ',num2str(ItCt),'...\n\n'));

        ImgStack_Input=ImgStack_WS_BVMEdge;
        ImgStack_PxIntensity_Input=ImgStack_Invert;
        ImgStack_Gradient_Input=ImgStack_Gradient;
        ImgStack_FgdMarker_Input=ImgStack_FgdMarker;

        for k=1:NumImgSlices

                % === Skip Z slices w/ 2x consistent WS
                if ~ismember(k,WS_It_ZFlagID)
                    WS_Prop_FlagID_Ct_Mtx(k,ItCt)=WS_Prop_FlagID_Ct_Mtx(k,ItCt-1);
                    continue;
                end
                
                % === Distinguish regions (potential cells) using regionprops
                [Img_WS_Lb,~]=bwlabel(ImgStack_Input(:,:,k));
                Img_PxIntensity=ImgStack_PxIntensity_Input(:,:,k);

                % Note: do regionprops after converting img to RGB; 
                % 1st dim=x, 2nd dim=y
                WS_regionprops=...
                    regionprops(Img_WS_Lb,Img_PxIntensity,...
                    'centroid','PixelIdxList','solidity','Perimeter','Area','ConvexHull','PixelValues');

                WS_Centroid_XY_List=cat(1,WS_regionprops.Centroid);
                WS_Centroid_LinIdx_List=sub2ind(size(ImgStack_Input(:,:,k)), round(WS_Centroid_XY_List(:,2)), round(WS_Centroid_XY_List(:,1)));

                WS_Solid_List=cat(1,WS_regionprops.Solidity);
                WS_Perim_List=cat(1,WS_regionprops.Perimeter);
                WS_Area_List=cat(1,WS_regionprops.Area);

                % Intensity: quantile difference
                WS_IntensityQtl_List=[];
                for j=1:size(WS_regionprops)
                    WS_Intensity_Quantile=...
                        quantile(WS_regionprops(j).PixelValues,[0.01,0.99]);     
                    WS_IntensityQtl_List=...
                        vertcat(WS_IntensityQtl_List,WS_Intensity_Quantile(1,2)-WS_Intensity_Quantile(1,1));
                end

                % === Assign regions with Candidate Cell ID 
                % Matrix has same format as Maxima_HullID, below       
                % Centroid_HullID=repmat((1:1:size(WS_regionprops,1))',1,2);

                % === Identify watershed regions disqualified as cells by bone, vessel and img masks
                WS_BoneDQ_FlagID=zeros(0,1);
                WS_VesselDQ_FlagID=zeros(0,1);
                WS_ImgMaskDQ_FlagID=zeros(0,1);
                WS_BCMaskDQ_FlagID=zeros(0,1);

                % Bone
                BoneMask=BoneMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
                WS_BoneMask_Flag=zeros(size(WS_regionprops,1),1);           
                for regionct=1:size(WS_regionprops)
                    Bone_OccpPixelCt=sum(BoneMask(WS_regionprops(regionct).PixelIdxList),1);
                    if Bone_OccpPixelCt/WS_Area_List(regionct)>0
                        WS_BoneMask_Flag(regionct)=1;                    
                    end
                end 
                WS_BoneDQ_FlagID=find(WS_BoneMask_Flag);

                % Vessel
                VesselMask=VesselMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
                WS_VesselMask_Flag=zeros(size(WS_regionprops,1),1);           
                for regionct=1:size(WS_regionprops)
                    Vessel_OccpPixelCt=sum(VesselMask(WS_regionprops(regionct).PixelIdxList),1);
                    if Vessel_OccpPixelCt/WS_Area_List(regionct)>0
                        WS_VesselMask_Flag(regionct)=1;                    
                    end
                end
                WS_VesselDQ_FlagID=find(WS_VesselMask_Flag);
                WS_VesselDQ_FlagID(ismember(WS_VesselDQ_FlagID,WS_BoneDQ_FlagID))=[]; % Assign to bone bone-vessel overlapping regions

                % BCMask
                if BCMask_Query=='y'
                    BCMask=BCMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
                    WS_BCMask_Flag=zeros(size(WS_regionprops,1),1);           
                    for regionct=1:size(WS_regionprops)
                        BC_OccpPixelCt=sum(~BCMask(WS_regionprops(regionct).PixelIdxList),1);
                        if BC_OccpPixelCt/WS_Area_List(regionct)>0
                            WS_BCMask_Flag(regionct)=1;                    
                        end
                    end
                    WS_BCMaskDQ_FlagID=find(WS_BCMask_Flag);
                    WS_BCMaskDQ_FlagID(ismember(WS_BCMaskDQ_FlagID,WS_BoneDQ_FlagID))=[]; % If conflict, assign to bone
                    WS_BCMaskDQ_FlagID(ismember(WS_BCMaskDQ_FlagID,WS_VesselDQ_FlagID))=[]; % If conflict, assign to vessels
                else
                    WS_BCMaskDQ_FlagID=[];
                end

                % ImgMask
                ImgMask=ImgMaskStack(:,:,k);    % = 1 matrix if no actual mask supplied
                WS_ImgMask_Flag=zeros(size(WS_regionprops,1),1);           
                for regionct=1:size(WS_regionprops)
                    ImgMask_OccpPixelCt=sum(~ImgMask(WS_regionprops(regionct).PixelIdxList),1);
                    if ImgMask_OccpPixelCt/WS_Area_List(regionct)>0
                        WS_ImgMask_Flag(regionct)=1;                    
                    end
                end
                WS_ImgMaskDQ_FlagID=find(WS_ImgMask_Flag);

                % === Identify regions disqualified as cells by regionprops

                % RULE 1. ID blobs with Area > Pix_Area_max or < Pix_Area_min
                if CellCt3D_Query=='y'
                    WS_Area_Flag=bsxfun(@gt,WS_Area_List,Px_Area_max)+bsxfun(@lt,WS_Area_List,Px_Area_min);
                    WS_Area_FlagID=find(WS_Area_Flag);        
                else
                    WS_Area_Flag=bsxfun(@gt,WS_Area_List,Px_Area_max)+bsxfun(@lt,WS_Area_List,Px_Area_min);
                    WS_Area_FlagID=find(WS_Area_Flag);        
                end
               WS_Area_FlagID(ismember(WS_Area_FlagID,...
                    unique(vertcat(WS_BoneDQ_FlagID,WS_VesselDQ_FlagID,WS_BCMaskDQ_FlagID,WS_ImgMaskDQ_FlagID))))=[];
%                     fprintf(strcat('... Slice',num2str(k),': ',num2str(numel(WS_Area_FlagID)),'/',num2str(size(WS_regionprops,1)),' WS areas to be removed w/ area rule...\n'));

                % RULE 2: Circularity (ImageJ definition; range 0-1)
% % %                 if CellCt3D_Query=='y'
% % %                     WS_Circ_FlagID=[];
% % %                 else
% % %                     WS_Circ_List=(4*pi.*WS_Area_List)./(WS_Perim_List.^2);
% % %                     WS_Circ_Flag=bsxfun(@lt,WS_Circ_List,0.5);
% % %                     WS_Circ_FlagID=find(WS_Circ_Flag);
% % %                 end
                WS_Circ_FlagID=[];

                WS_Circ_FlagID(ismember(WS_Circ_FlagID,...
                    unique(vertcat(WS_BoneDQ_FlagID,WS_VesselDQ_FlagID,WS_BCMaskDQ_FlagID,WS_ImgMaskDQ_FlagID))))=[];
%                     fprintf(strcat('... Slice',num2str(k),': ',num2str(numel(WS_Circ_FlagID)),'/',num2str(size(WS_regionprops,1)),' WS areas to be removed w/ circularity rule...\n'));

                % RULE 3: Solidity (Area/CvxHullArea; range 0-1)
% % %                 if CellCt3D_Query=='y'
% % %                     WS_Solid_Flag=bsxfun(@lt,WS_Solid_List,0.5);
% % %                     WS_Solid_FlagID=[];
% % %                 else
% % %                     WS_Solid_Flag=bsxfun(@lt,WS_Solid_List,0.5);  % Default 0.7
% % %                     WS_Solid_FlagID=find(WS_Solid_Flag);
% % %                 end
                WS_Solid_FlagID=[];
                WS_Solid_FlagID(ismember(WS_Solid_FlagID,...
                    unique(vertcat(WS_BoneDQ_FlagID,WS_VesselDQ_FlagID,WS_BCMaskDQ_FlagID,WS_ImgMaskDQ_FlagID))))=[];
%                     fprintf(strcat('... Slice',num2str(k),': ',num2str(numel(WS_Solid_FlagID)),'/',num2str(size(WS_regionprops,1)),' WS areas to be removed w/ solidity rule...\n'));

                % RULE 4: Brightness
                if CellCt3D_Query=='y'
                    WS_IntensityQtl_Flag=bsxfun(@lt,WS_IntensityQtl_List,BC_BkgdNoiseRange); 
                    WS_IntensityQtl_FlagID=find(WS_IntensityQtl_Flag);
                else
                    WS_IntensityQtl_Flag=bsxfun(@lt,WS_IntensityQtl_List,BC_BkgdNoiseRange);
                    WS_IntensityQtl_FlagID=find(WS_IntensityQtl_Flag);
                end
               WS_IntensityQtl_FlagID(ismember(WS_IntensityQtl_FlagID,...
                    unique(vertcat(WS_BoneDQ_FlagID,WS_VesselDQ_FlagID,WS_BCMaskDQ_FlagID,WS_ImgMaskDQ_FlagID))))=[];
%                     fprintf(strcat('... Slice',num2str(k),': ',num2str(numel(WS_IntensityQtl_FlagID)),'/',num2str(size(WS_regionprops,1)),' WS areas to be removed w/ Intensity rule...\n'));

                % Combine all flags
                WS_Prop_FlagID=zeros(0,1);
                WS_Prop_FlagID=union(union(union(WS_Area_FlagID,WS_Circ_FlagID),WS_Solid_FlagID),WS_IntensityQtl_FlagID);
                fprintf(strcat('... Slice',num2str(k),': ',num2str(numel(WS_Prop_FlagID)),'/',num2str(size(WS_regionprops,1)),' WS areas to be removed...\n'));
                WS_Prop_FlagID_Ct_Mtx(k,ItCt)=numel(WS_Prop_FlagID);
                
                % === Erase "unqualified" Watershed regions        
                % Combine all disqualified FlagID; erase from lists
                WS_DQ_FlagID=vertcat(reshape(WS_BoneDQ_FlagID,numel(WS_BoneDQ_FlagID),1),...
                    reshape(WS_VesselDQ_FlagID,numel(WS_VesselDQ_FlagID),1),...
                    reshape(WS_BCMaskDQ_FlagID,numel(WS_BCMaskDQ_FlagID),1),...
                    reshape(WS_ImgMaskDQ_FlagID,numel(WS_ImgMaskDQ_FlagID),1),...
                    reshape(WS_Prop_FlagID,numel(WS_Prop_FlagID),1));
                WS_DQ_FlagID=unique(WS_DQ_FlagID);
                % WS_DQ_FlagID=union(union(union(WS_BoneDQ_FlagID,WS_VesselDQ_FlagID),WS_ImgMaskDQ_FlagID),WS_Prop_FlagID);

                WS_Cell_FlagID=(1:1:size(WS_regionprops,1))';
                WS_Cell_FlagID(ismember(WS_Cell_FlagID,WS_DQ_FlagID))=[];

                WS_Cell_Centroid_XY_List=WS_Centroid_XY_List(WS_Cell_FlagID,:);
                WS_Cell_Centroid_LinIdx_List=WS_Centroid_LinIdx_List(WS_Cell_FlagID,:);
                WS_Cell_Area_List=WS_Area_List(WS_Cell_FlagID,:);

                % === Output
                ImgStack_WS_Lb(:,:,k)=Img_WS_Lb;
                WS_regionprops_Cell{k,1}=WS_regionprops;

                WS_BoneDQ_FlagID_Cell{k,1}=WS_BoneDQ_FlagID;
                WS_VesselDQ_FlagID_Cell{k,1}=WS_VesselDQ_FlagID;
                WS_BCMaskDQ_FlagID_Cell{k,1}=WS_BCMaskDQ_FlagID;
                WS_ImgMaskDQ_FlagID_Cell{k,1}=WS_ImgMaskDQ_FlagID;
                WS_PropDQ_FlagID_Cell{k,1}=WS_Prop_FlagID;
                WS_DQ_FlagID_Cell{k,1}=WS_DQ_FlagID;

                WS_Cell_FlagID_Cell{k,1}=WS_Cell_FlagID;

                WS_Cell_Centroid_XY_List_Cell{k,1}=WS_Cell_Centroid_XY_List;  
                WS_Cell_Centroid_LinIdx_List_Cell{k,1}=WS_Cell_Centroid_LinIdx_List;  
                WS_Cell_Area_List_Cell{k,1}=WS_Cell_Area_List;  
                
                % === Clean up WS Foreground Markers & Gradient: 
                % Dim regions designated as DQ regions for this step
                % * Set foreground markers to 0 for subset of DQ regions
                % * Set gradient of pixels within DQ regions to 0
   
                MCWS_FlagID_PixelIdxList=cat(1,WS_regionprops(WS_IntensityQtl_FlagID).PixelIdxList);

                Img_Gradient_PostWS=ImgStack_Gradient_Input(:,:,k);        
                Img_Gradient_PostWS(MCWS_FlagID_PixelIdxList)=0; % Set gradient = 0

                Img_FgdMarker_PostWS=ImgStack_FgdMarker_Input(:,:,k);        
                Img_FgdMarker_PostWS(MCWS_FlagID_PixelIdxList)=0; % Remove foreground markers

                [PostWS_Lb,~]=bwlabel(Img_Gradient_PostWS); % Also remove partially removed foreground markers
                PostWS_regionprops=regionprops(PostWS_Lb,'Area','PixelIdxList');
                PostWS_Area_List=cat(1,PostWS_regionprops.Area);
                PostWS_Area_Flag=bsxfun(@lt,PostWS_Area_List,FgdMarker_Area_List_Min);
                PostWS_Area_FlagID=find(PostWS_Area_Flag);   
                FgdMarkerPostWS_Area_FlagID_PixelIdxList=...
                    cat(1,PostWS_regionprops(PostWS_Area_FlagID).PixelIdxList);
                 Img_Gradient_PostWS(FgdMarkerPostWS_Area_FlagID_PixelIdxList)=0;

                % Prepare for next MCWS Iteration
                ImgStack_Gradient_Input(:,:,k)=Img_Gradient_PostWS;
                ImgStack_FgdMarker(:,:,k)=Img_FgdMarker_PostWS;

        end
        
    % Create matrix designating which Z slice requires more iteration
    % Stop iteration for each Z if its WS_Prop_FlagID_Ct is consistent for 2 rounds
    if ItCt>2
        WS_It_ZFlagID_Old=WS_It_ZFlagID;
        WS_It_ZFlagID=find(logical(...
            abs(WS_Prop_FlagID_Ct_Mtx(:,ItCt)-WS_Prop_FlagID_Ct_Mtx(:,ItCt-1))+...
            abs(WS_Prop_FlagID_Ct_Mtx(:,ItCt)-WS_Prop_FlagID_Ct_Mtx(:,ItCt-2))...
            ));
        WS_It_ZFlagID_EqualCheck=WS_It_ZFlagID_EqualCheck+isequal(WS_It_ZFlagID_Old,WS_It_ZFlagID);
        if WS_It_ZFlagID_EqualCheck>5
            fprintf(strcat('\nTerminate Iteration: unable to stabilize WS_Prop_FlagID_Ct. Iteration ',num2str(ItCt),'...\n\n'));
            break;
        end
        if isempty(WS_It_ZFlagID) % All done
            ItCt=-9999;
        else
            ItCt=ItCt+1;
        end

    else
        ItCt=ItCt+1;    % First 2 iterations
    end

end

WS_Prop_FlagID_Ct_Mtx(:,find(isnan(WS_Prop_FlagID_Ct_Mtx(1,:))))=[];

clearvars ImgStack_OBRCBR;
clearvars Img_* -except Img_Height Img_Width;
clearvars FgdMarker* -except  FgdMarker_Area_List_Min;
clearvars ImgStack*Gradient;
clearvars ImgStack*Marker -except ImgStack_FgdMarker ImgStack_BkgdMarker;

clearvars ImgStack_WS* -except ImgStack_WS_Lb ImgStack_WS_BVMEdge;
clearvars ImgStack_Gradient*;
clearvars WS_* -except WS_*_Cell WS_Prop_FlagID_Ct_Mtx;
clearvars PostWS* MCWS*;
clearvars Img_* -except Img_Height Img_Width;
clearvars MaskStackEdge;
clearvars WS_It_ZFlagID*;

clearvars BC_Bkgd*;
clearvars *_Input;


%% Save Images

fprintf('Saving Watershed-related Images...\n');

for k=1:NumImgSlices
    
    % Pre Watershed Stack (Image Stack used for Watershed) (type 'double')
    % imwrite(double(ImgStack_Invert(:,:,k)),strcat(SavePreWSFilePath,'PreWS -',num2str(k,'%04.0f'),'.tif'));
    
    % Watershed Stack (Image Stack used for Watershed) (type 'double')
    % imwrite(double(ImgStack_WS_BVMEdge(:,:,k)),strcat(SaveWSLineFilePath,'WSLine -',num2str(k,'%04.0f'),'.tif'));
    
    % Foreground Marker
    % imwrite(double(ImgStack_FgdMarker(:,:,k)),strcat(SaveWSFgdMFilePath,'WSFgd -',num2str(k,'%04.0f'),'.tif'));

    % Background Marker
    % imwrite(double(ImgStack_BkgdMarker(:,:,k)),strcat(SaveWSBgdMFilePath,'WSBgd -',num2str(k,'%04.0f'),'.tif'));
    
end

clearvars ImgStack_WS_BVMEdge;
clearvars ImgStack_FgdMarker ImgStack_BkgdMarker;
clearvars *_Input;


%% Prepare "Checkerboard" Grid for Bone and Vessel Localization
% Grid ensures non-WS cell areas (bone, vessel, Img masks) have proper
% boundaries; better area and density measurements

fprintf('"Checkerboard" Grid for area and density measurements...\n');

% Input
WS_regionprops_Cell_Input=WS_regionprops_Cell;
WS_DQ_FlagID_Cell_Input=WS_DQ_FlagID_Cell;

% Output
Gd_Cell_FlagID_Cell=cell(NumImgSlices,1);
Gd_BoneDQ_FlagID_Cell=cell(NumImgSlices,1);
Gd_VesselDQ_FlagID_Cell=cell(NumImgSlices,1);
Gd_BCMaskDQ_FlagID_Cell=cell(NumImgSlices,1);
Gd_ImgMaskDQ_FlagID_Cell=cell(NumImgSlices,1);
    
% === Create Grid, label each checker and analyze properties
% Grid is the same across all Z Slice
Img_Gd=ones(Img_Height,Img_Width,'logical');
Img_Gd(Kernel_px_min:Kernel_px_min:Kernel_px_min*floor(Img_Height/Kernel_px_min),:)=0;
Img_Gd(:,Kernel_px_min:Kernel_px_min:Kernel_px_min*floor(Img_Width/Kernel_px_min))=0;

[Img_Gd_Lb, Img_Gd_NumRegions]=bwlabel(Img_Gd);   % Label checkers

% Obtain Grid checker properties
Gd_regionprops=regionprops(Img_Gd_Lb,'PixelIdxList','Area');
Gd_Area_List=cat(1,Gd_regionprops.Area);

% === Identify grid regions occupied by cell, bone, vessel and img masks

for k=1:size(WS_regionprops_Cell_Input,1)    

        Gd_WSCell_FlagID=zeros(0,1);
        Gd_BoneDQ_FlagID=zeros(0,1);
        Gd_VesselDQ_FlagID=zeros(0,1);
        Gd_BCMaskDQ_FlagID=zeros(0,1);
        Gd_ImgMaskDQ_FlagID=zeros(0,1);

        % Cell
        % Create Cell Mask by removing WS regions disqualified by bone, vessel, ImgMask
        WS_regionprops=WS_regionprops_Cell_Input{k,1};  
        WS_DQ_FlagID=WS_DQ_FlagID_Cell_Input{k,1};
        
        WS_regionprops(WS_DQ_FlagID)=[];
        WS_PixelIdxList_List=cat(1,WS_regionprops(:).PixelIdxList);
        CellMask=zeros(Img_Height,Img_Width,'logical');
        CellMask(WS_PixelIdxList_List)=1;    % = 0 matrix if no cells present

        Gd_CellDQ_Flag=zeros(Img_Gd_NumRegions,1);           
        for gridct=1:Img_Gd_NumRegions
            Cell_OccpPixelCt=sum(CellMask(Gd_regionprops(gridct).PixelIdxList),1);
            if Cell_OccpPixelCt/Gd_Area_List(gridct)>0  
                Gd_CellDQ_Flag(gridct)=1;                    
            end
        end 
        Gd_CellDQ_FlagID=find(Gd_CellDQ_Flag);

        % Bone
        BoneMask=BoneMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
        Gd_BoneDQ_Flag=zeros(Img_Gd_NumRegions,1);           
        for gridct=1:Img_Gd_NumRegions
            Bone_OccpPixelCt=sum(BoneMask(Gd_regionprops(gridct).PixelIdxList),1);
            if Bone_OccpPixelCt/Gd_Area_List(gridct)>0.5
                Gd_BoneDQ_Flag(gridct)=1;                    
            end
        end 
        Gd_BoneDQ_FlagID=find(Gd_BoneDQ_Flag);
        Gd_BoneDQ_FlagID(ismember(Gd_BoneDQ_FlagID,Gd_CellDQ_FlagID))=[]; % exclude cell grids

        % Vessel
        VesselMask=VesselMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
        Gd_VesselDQ_Flag=zeros(Img_Gd_NumRegions,1);           
        for gridct=1:Img_Gd_NumRegions
            Vessel_OccpPixelCt=sum(VesselMask(Gd_regionprops(gridct).PixelIdxList),1);
            if Vessel_OccpPixelCt/Gd_Area_List(gridct)>0.5
                Gd_VesselDQ_Flag(gridct)=1;                    
            end
        end
        Gd_VesselDQ_FlagID=find(Gd_VesselDQ_Flag);
        Gd_VesselDQ_FlagID(ismember(Gd_VesselDQ_FlagID,Gd_CellDQ_FlagID))=[]; % exclude cell grids
        Gd_VesselDQ_FlagID(ismember(Gd_VesselDQ_FlagID,Gd_BoneDQ_FlagID))=[]; % Assign to bone bone-vessel overlapping regions
        
        % BCMask
        BCMask=BCMaskStack(:,:,k);  % = 0 matrix if no actual mask supplied
        Gd_BCMaskDQ_Flag=zeros(Img_Gd_NumRegions,1);           
        for gridct=1:Img_Gd_NumRegions
            BCMask_OccpPixelCt=sum(~BCMask(Gd_regionprops(gridct).PixelIdxList),1);
            if BCMask_OccpPixelCt/Gd_Area_List(gridct)>0.5
                Gd_BCMaskDQ_Flag(gridct)=1;                    
            end
        end
        Gd_BCMaskDQ_FlagID=find(Gd_BCMaskDQ_Flag);
        Gd_BCMaskDQ_FlagID(ismember(Gd_BCMaskDQ_FlagID,Gd_CellDQ_FlagID))=[]; % exclude cell grids
        Gd_BCMaskDQ_FlagID(ismember(Gd_BCMaskDQ_FlagID,Gd_BoneDQ_FlagID))=[]; % If conflict, assign to bone
        Gd_BCMaskDQ_FlagID(ismember(Gd_BCMaskDQ_FlagID,Gd_VesselDQ_FlagID))=[]; % If conflict, assign to vessel

        % ImgMask
        ImgMask=ImgMaskStack(:,:,k);    % = 1 matrix if no actual mask supplied
        Gd_ImgMaskDQ_Flag=zeros(Img_Gd_NumRegions,1);           
        for gridct=1:Img_Gd_NumRegions
            ImgMask_OccpPixelCt=sum(~ImgMask(Gd_regionprops(gridct).PixelIdxList),1);
            if ImgMask_OccpPixelCt/Gd_Area_List(gridct)>0.5
                Gd_ImgMaskDQ_Flag(gridct)=1;                    
            end
        end
        Gd_ImgMaskDQ_FlagID=find(Gd_ImgMaskDQ_Flag);

        % === Output
        Gd_Cell_FlagID_Cell{k,1}=Gd_CellDQ_FlagID;
        Gd_BoneDQ_FlagID_Cell{k,1}=Gd_BoneDQ_FlagID;
        Gd_VesselDQ_FlagID_Cell{k,1}=Gd_VesselDQ_FlagID;
        Gd_BCMaskDQ_FlagID_Cell{k,1}=Gd_BCMaskDQ_FlagID;
        Gd_ImgMaskDQ_FlagID_Cell{k,1}=Gd_ImgMaskDQ_FlagID;
        % NOTE: Img_Gd_Lb is also output
        
end

clearvars WS_PixelIdxList_List;
clearvars BoneMaskStack VesselMaskStack BCMaskStack;
clearvars Img_Gd Img_Gd_NumRegions;
clearvars Gd_Area_List Gd_regionprops;
clearvars Gd*FlagID;
clearvars *Occp*;
clearvars CellMask BoneMask* VesselMask* BCMask* ImgMask*;
clearvars ImgMaskTranslate;
clearvars *_Input;


%% Candidate Cell Identification using "Find Maxima" Function in ImageJ
% "Find Maxima" locates center of candidate cells, no boundaries
%*****

fprintf('Candidate Cell Identification by ImageJ "Find Maxima"...\n');

ImgStack_Input=ImgStack_Invert;

% Cell for holding Maxima_XY_List for all slices
Maxima_XY_List_Cell=cell(size(ImgStack_Input,3),1);

if FindMaxima_Query=='y'
    
    % Start Miji
    javaaddpath '/Applications/MATLAB_R2017b.app/java/mij.jar';
    javaaddpath '/Applications/MATLAB_R2017b.app/java/ij-1.51s.jar';

    MIJ.start;

    for k = 1:size(ImgStack_Input,3)
    % for k = 1:1

           Img=single(ImgStack_Input(:,:,k));

           % === Scale for ImageJ 
           PreFMScaleFactor=(2^31-1)/max(Img(:));
           Img_PreFMScale=single(Img.*PreFMScaleFactor);

            % === Find Maxima @ ImageJ
            MIJ.run('Close All');
            MIJ.createImage(Img_PreFMScale);
            MIJ.run("32-bit"); 

            % FindMaxima_commandstr=sprintf(strcat('\''','Find Maxima...','\''',',','\''','noise=',num2str(10),' output=List exclude','\'''));
            % MIJ.run(FindMaxima_commandstr);
            % MIJ.run('Find Maxima...', 'noise=10 output=[Segmented Particles] exclude');
            % MIJ.run('Find Maxima...', 'noise=10 output=[Single Points] exclude');
            MIJ.run('Find Maxima...', 'noise=0 output=List exclude');  
            pause(1);

            Maxima_XY_List=MIJ.getResultsTable;
            MIJ.run('Clear Results');

            if isempty(Maxima_XY_List)
                Maxima_XY_List=[-9999 -9999];
            end

            Maxima_XY_List_Cell{k,1}=Maxima_XY_List;

    end

    MIJ.run('Close All');
    MIJ.exit;

elseif FindMaxima_Query=='n'
    
    for k = 1:size(ImgStack_Input,3)
        Maxima_XY_List_Cell{k,1}=[-9999 -9999];
    end
    
end

clearvars ImgStack_Invert;
clearvars *_Input;


%% Cell identification from Region IDs from "Find Maxima" Function in ImageJ

% As "Find Maxima" locates center of cell (no boundaries),
% Use MATLAB watershed info for cell boundaries.
% Locate watershed region for each Maxima using inhull        
% Create matrix that records which convex hull (col) maxima (row) is found
    
fprintf('Fitting "Find Maxima" pts into Watershed region convex hulls ...\n');

Maxima_XY_List_Cell_Input=Maxima_XY_List_Cell;
WS_regionprops_Cell_Input=WS_regionprops_Cell;
WS_DQ_FlagID_Cell_Input=WS_DQ_FlagID_Cell;

% Prepare Cells for output:
% 1) Matrix matching each FindMaxima Pt (Col1) to Watershed regions (Col2)
% 2) List of FindMaxima pt (in linear index) occupied Watershed regions identified as cells, 
% after regionprops and bone/vessel/ImgMask disqualification
% 3) Area of FindMaxima pt occupied Watershed regions identified as cells, 
% after regionprops and bone/vessel/ImgMask disqualification
% 4) List of FindMaxima pt (in linear index) without watershed region assigned
Maxima_Cell_WSHullID_Cell=cell(NumImgSlices,1);
Maxima_Cell_XY_List_Cell=cell(NumImgSlices,1);
Maxima_Cell_LinIdx_List_Cell=cell(NumImgSlices,1);
Maxima_Cell_Area_List_Cell=cell(NumImgSlices,1);
Maxima_NoWSHull_LinIdx_List_Cell=cell(NumImgSlices,1);

if FindMaxima_Query=='y'
    for k=1:NumImgSlices

        Maxima_XY_List=Maxima_XY_List_Cell_Input{k,1};
        WS_regionprops=WS_regionprops_Cell_Input{k,1};
        WS_DQ_FlagID=WS_DQ_FlagID_Cell_Input{k,1};

        % === Change MaximaXY_List to Linear Indices
        % X=col, left to right; Y=row, top to bottom
        if Maxima_XY_List(1,1)==-9999
            Maxima_LinIdx_List=-9999;
        else
            Maxima_LinIdx_List=sub2ind(size(Img), Maxima_XY_List(:,2), Maxima_XY_List(:,1));
        end

        % === Match each Maxima to Watershed regions using inhull   
        % "Maxima_WSHullID" has 2 columns:
        % Col1: index of Maxima pt (1:1:Number of Maxima Pts)
        % Col2: WS Region (as Convex Hull) number in which Maxima pt in Col1 is located        
        Maxima_MultiWSHull_Idx=zeros(0,1);
        Maxima_NoWSHull_Idx=zeros(0,1);

        if Maxima_LinIdx_List(1)==-9999
            Maxima_WSHullID=[-9999 -9999];
        else
            inhull_idx=zeros(size(Maxima_LinIdx_List,1),size(WS_regionprops,1));

            % ConvexHull indices also given in col1=X col2=Y, 
            % same as Maxima_XY_List
            for hullct=1:size(WS_regionprops,1)
                inhull_idx(:,hullct)=inhull(Maxima_XY_List,WS_regionprops(hullct).ConvexHull,[],Kernel_px_min);
            end        

            % Determine which Maxima is missing / have multiple blob assignment        
            % For Maxima with multiple hull assignment, assign to WS
            % region with shortest distance centroid 
            if Maxima_LinIdx_List(1)~=-9999
                Hull_OccCt=sum(inhull_idx,2);   % Col vector; each row for one MaximaXY pt
                Maxima_MultiWSHull_Idx=find(Hull_OccCt>1);
                Maxima_NoWSHull_Idx=find(Hull_OccCt==0);
            end

            if ~isempty(Maxima_MultiWSHull_Idx)
                for i=1:size(Maxima_MultiWSHull_Idx,1)
                % for i=1:1
                    Hull=find(inhull_idx(Maxima_MultiWSHull_Idx(i),:))';
                    Hull_CentroidXY=cat(1,WS_regionprops(Hull).Centroid);    % Each row=centroid XY for each hull
                    [~,HullMinDistIdx]=min((Hull_CentroidXY(:,1)-Maxima_XY_List(Maxima_MultiWSHull_Idx(i),1)).^2+(Hull_CentroidXY(:,2)-Maxima_XY_List(Maxima_MultiWSHull_Idx(i),2)).^2);
                    Hull(HullMinDistIdx)=[];    % Remove minimal distance hull (to be kept) from Hull List
                    inhull_idx(Maxima_MultiWSHull_Idx(i),Hull)=0;   % Set non-min centroid-to-MaximaXY distance hull to 0
                end
            end

            [inhull_idx_NZrow,inhull_idx_NZcol]=find(inhull_idx);   % Row=Maxima Pt; Col=Corresponding Hull (WS region)
            Maxima_WSHullID=sortrows([inhull_idx_NZrow,inhull_idx_NZcol],1);  

            % Determine which WS region (Hull) have multiple maxima assignment        
            % For WS regions with multiple maxima assignment, assign to 
            % Maxima with the shortest distance to WS region's centroid          
            if Maxima_WSHullID(1)~=-9999
                HistCtEdge=1:1:max(Maxima_WSHullID(:,2))+1; 
                [WSHull_MaximaHistCt,~] = histcounts(Maxima_WSHullID(:,2),HistCtEdge);
                WSHull_MultiMax_Idx=find(WSHull_MaximaHistCt>1);    % WS regions (Hulls) occupied by more than one maxima
                if ~isempty(WSHull_MultiMax_Idx)
                    for j=1:size(WSHull_MultiMax_Idx,1)
                        Maxima_Idx=Maxima_WSHullID(Maxima_WSHullID(:,2)==WSHull_MultiMax_Idx(j),1);
                        Maxima_XY=Maxima_XY_List(Maxima_Idx,:); % Each row=Maxima XY sharing WS region (Hull)
                        WSHull_CentroidXY=WS_regionprops(Maxima_MultiWSHull_Idx(j)).Centroid;
                        [~,MaximaMinDistIdx]=min((Maxima_XY(:,1)-WSHull_CentroidXY(:,1)).^2+(Maxima_XY(:,2)-WSHull_CentroidXY(:,2)).^2);
                        Maxima_Idx(MaximaMinDistIdx)=[];    % Remove Maxima further away from WS region centroid
                        Maxima_WSHullID(ismember(Maxima_WSHullID(:,1),Maxima_Idx),:)=[];
                    end
                end
            end

        end

        % === Erase "unqualified" Watershed regions        
        if Maxima_LinIdx_List(1)==-9999
            Maxima_Cell_WSHullID=[-9999 -9999];
            Maxima_Cell_XY_List=[-9999 -9999];
            Maxima_Cell_LinIdx_List=-9999;
            Maxima_Cell_Area_List=-9999;
        else
            Maxima_Cell_WSHullID=Maxima_WSHullID;
            Maxima_Cell_WSHullID(ismember(Maxima_Cell_WSHullID(:,2),WS_DQ_FlagID),:)=[];
            Maxima_Cell_XY_List=Maxima_XY_List(Maxima_Cell_WSHullID(:,1),:);  % XY coordinates of Maxima pt of Cells  
            Maxima_Cell_LinIdx_List=Maxima_LinIdx_List(Maxima_Cell_WSHullID(:,1));  % Linear idx of Maxima pt of Cells          
            WS_Area_List=cat(1,WS_regionprops.Area);    % Area of Maxima Cells
            Maxima_Cell_Area_List=WS_Area_List(Maxima_Cell_WSHullID(:,2));
            if isempty(Maxima_Cell_LinIdx_List)
                Maxima_Cell_WSHullID=[-9999 -9999];
                Maxima_Cell_LinIdx_List=-9999;  % Maxima coordiate as 
                Maxima_Area_LinIdx_List=-9999;              
            end       
        end

        % === Output
        % List of Cells identified by FindMaxima + Watershed, after region and mask
        % disqualification       
        Maxima_Cell_WSHullID_Cell{k,1}=Maxima_Cell_WSHullID;
        Maxima_Cell_XY_List_Cell{k,1}=Maxima_Cell_XY_List;
        Maxima_Cell_LinIdx_List_Cell{k,1}=Maxima_Cell_LinIdx_List;
        Maxima_Cell_Area_List_Cell{k,1}=Maxima_Cell_Area_List;      
        if Maxima_LinIdx_List(1)~=-9999
            Maxima_NoWSHull_LinIdx_List_Cell{k,1}=Maxima_LinIdx_List(Maxima_NoWSHull_Idx');
        else
            Maxima_NoWSHull_LinIdx_List_Cell{k,1}=-9999;
        end

    end
elseif FindMaxima_Query=='n'
    Maxima_Cell_WSHullID_Cell(:)={[-9999 -9999]};
    Maxima_Cell_XY_List_Cell(:)={[-9999 -9999]};
    Maxima_Cell_LinIdx_List_Cell(:)={-9999};
    Maxima_Cell_Area_List_Cell(:)={-9999};   
    Maxima_NoWSHull_LinIdx_List_Cell(:)={-9999};
end

clearvars WS_regionprops;
clearvars Maxima_XY_List_Cell WS_regionprops_Cell WS_DQ_FlagID_Cell;
clearvars *_Input;


%%  Brightness Contrast Adjustment for local image saving. 
% Do not affect calculations. Do Not Remove.

ImgStack_Temp=reshape(ImgStack,numel(ImgStack),1,1);
ImgStack_Temp(ImgStack_Temp<0.2)=[];
ImgStack_Temp(ImgStack_Temp>1-eps)=[];
ImgStack_BCAdjMean=mean(ImgStack_Temp(:));

clearvars ImgStack_Temp;


%% Colocalization check: Cells as coordinate pts on Images

fprintf('Colocalization check (ImgStack, OrgImgStack (Enhanced Contrast)): Cells as coordinates ...\n');

% INPUT: Linear index of cell point coordinate
WS_Cell_Centroid_LinIdx_List_Cell_Input=WS_Cell_Centroid_LinIdx_List_Cell;       % METHOD1: Watershed Centroid Method
Maxima_Cell_LinIdx_List_Cell_Input=Maxima_Cell_LinIdx_List_Cell;     % METHOD2: Find Maxima Method
        
% === Prepare stack with cell pts = 1, rest = 0

Dilate_SE = strel('disk',2);    % Dilation for easier to read pts

% Prepare output cells
% 3 Columns for 2 METHODs
% 1) Colocalization: WSCentroid(R) on ImgStack
% 2) Colocalization: Maxima(G) on OrgImgStack
% 3) Colocalization: WSCentroid(R), Maxima(G) on OrgImgStack
ImgRGB_WSCentroid_Maxima_Colocal_Cell=cell(size(WS_Cell_Centroid_LinIdx_List_Cell_Input,1),2);
OrgImgRGB_WSCentroid_Maxima_Colocal_Cell=cell(size(WS_Cell_Centroid_LinIdx_List_Cell_Input,1),2);

for METHOD=1:2  
                 
        % === Create colocal images
        for k=1:size(WS_Cell_Centroid_LinIdx_List_Cell_Input,1)
            
                % === Choose Method
                if METHOD==1
                    CellPtLinList=WS_Cell_Centroid_LinIdx_List_Cell_Input{k,1};
                elseif METHOD==2
                    CellPtLinList=Maxima_Cell_LinIdx_List_Cell_Input{k,1};
                end
                      
                % Create slice with cell pts = 255, rest 0
                CellPtSlice=zeros(Img_Height,Img_Width,'logical');               
                if (~isempty(CellPtLinList))    % Watershed regions present
                    if (CellPtLinList(1)~=-9999)  % -9999 only if no Find Maxima cell pts
                        CellPtSlice(CellPtLinList)=1;
                    end                
                    CellPtSlice=imdilate(CellPtSlice,Dilate_SE);
                end

                % Overlay on ImgStack slice
                % Color assignment: R=METHOD1; G=METHOD2
                CellPtImgOverlay=max(CellPtSlice,ImgStack(:,:,k));
                ImgRGB_CellPt_Colocal=repmat(ImgStack(:,:,k).*1/ImgStack_BCAdjMean,1,1,3); 
                ImgRGB_CellPt_Colocal=ImgRGB_CellPt_Colocal.*1;  % 2x intensity for visualization
                ImgRGB_CellPt_Colocal(ImgRGB_CellPt_Colocal>1)=1;   % 2x intensity for visualization
                ImgRGB_CellPt_Colocal(:,:,METHOD)=CellPtImgOverlay;  
                ImgRGB_WSCentroid_Maxima_Colocal_Cell{k,METHOD}=ImgRGB_CellPt_Colocal;                
                if METHOD==2
                    ImgRGB_WSCentroid_Maxima_Colocal_Cell{k,3}=...  % R (METHOD1) + G (METHOD2) 
                        max(ImgRGB_WSCentroid_Maxima_Colocal_Cell{k,1},ImgRGB_WSCentroid_Maxima_Colocal_Cell{k,2});
                end

                % Overlay on OrgImgStack slice
                % Color assignment: R=METHOD1; G=METHOD2
                Im_Layer=imadjust(OrgImgStack(:,:,k),stretchlim(OrgImgStack(:,:,k)));
                CellPtOrgImgOverlay=max(CellPtSlice,Im_Layer);
                OrgImgRGB_CellPt_Colocal=repmat(Im_Layer,1,1,3); 
                OrgImgRGB_CellPt_Colocal=OrgImgRGB_CellPt_Colocal*1;  % 2x intensity for visualization
                OrgImgRGB_CellPt_Colocal(OrgImgRGB_CellPt_Colocal>1)=1;   % 2x intensity for visualization
                OrgImgRGB_CellPt_Colocal(:,:,METHOD)=CellPtOrgImgOverlay;  
                OrgImgRGB_WSCentroid_Maxima_Colocal_Cell{k,METHOD}=OrgImgRGB_CellPt_Colocal;                    
                 if METHOD==2
                    OrgImgRGB_WSCentroid_Maxima_Colocal_Cell{k,3}=...  % R (METHOD1) + G (METHOD2) 
                        max(OrgImgRGB_WSCentroid_Maxima_Colocal_Cell{k,1},OrgImgRGB_WSCentroid_Maxima_Colocal_Cell{k,2});
                end


        end

end

clearvars CellPtLinList;
clearvars CellPtSlice;
clearvars ImgRGB_CellPt_Colocal CellPtImgOverlay;
clearvars OrgImgRGB_CellPt_Colocal CellPtOrgImgOverlay;
clearvars Im_Layer;
clearvars *_Input;


%% Save Images
% Colocalization, cell coordinates w/ ImgStack & OrgImgStack
% RED: Cell Watershed Centroid; GREEN: Cell Maxima 

fprintf('Saving Cell Coordinate Colocalization Images...\n');

for k=1:NumImgSlices

    % Colocalization with ImgStack
    % imwrite(double(ImgRGB_WSCentroid_Maxima_Colocal_Cell{k,3}),strcat(SaveColocalCellPosImgFilePath,'ColocalCellPosImg -',num2str(k,'%04.0f'),'.tif'));

    % Colocalization, w/ OrgImgStack
    imwrite(double(OrgImgRGB_WSCentroid_Maxima_Colocal_Cell{k,3}),strcat(SaveColocalCellPosOrgImgFilePath,'ColocalCellPosOrgImg -',num2str(k,'%04.0f'),'.tif'));           

end

clearvars ImgRGB_WSCentroid_Maxima_Colocal_Cell OrgImgRGB_WSCentroid_Maxima_Colocal_Cell;


%% Prepare Color Map of Cells, Bone, Vessel, ImgMask, Regionprops DQed regions for colocalization check:
% Assign different color index to cells, bones, vessels etc, 
% then assign color to color index

fprintf('Colormaps of Cells, Bone, Vessel, ImgMask for colocalization check ...\n');

ImgStack_WS_Lb_Input=ImgStack_WS_Lb;    %  Label by Watershed regions
WS_BoneDQ_FlagID_Cell_Input=WS_BoneDQ_FlagID_Cell;
WS_VesselDQ_FlagID_Cell_Input=WS_VesselDQ_FlagID_Cell;
WS_BCMaskDQ_FlagID_Cell_Input=WS_BCMaskDQ_FlagID_Cell;
WS_ImgMaskDQ_FlagID_Cell_Input=WS_ImgMaskDQ_FlagID_Cell;
WS_Prop_FlagID_Cell_Input=WS_PropDQ_FlagID_Cell;

Maxima_Cell_WSHullID_Cell_Input=Maxima_Cell_WSHullID_Cell;

Img_Gd_Lb_Input=Img_Gd_Lb;
Gd_CellDQ_FlagID_Cell_Input=Gd_Cell_FlagID_Cell;        
Gd_BoneDQ_FlagID_Cell_Input=Gd_BoneDQ_FlagID_Cell;
Gd_VesselDQ_FlagID_Cell_Input=Gd_VesselDQ_FlagID_Cell;
Gd_BCMaskDQ_FlagID_Cell_Input=Gd_BCMaskDQ_FlagID_Cell;
Gd_ImgMaskDQ_FlagID_Cell_Input=Gd_ImgMaskDQ_FlagID_Cell;

% Output
ImgStack_WS_Lb_IdxMap=zeros(size(ImgStack),'uint8');
ImgStack_Gd_Lb_IdxMap=zeros(size(ImgStack),'uint8');
ImgStack_Gd_Lb_NoEdgeIdxMap=zeros(size(ImgStack),'uint8');

Img_WS_Lb_MapRGB_Cell=cell(size(ImgStack_WS_Lb_Input,3),1);
Img_Gd_Lb_MapRGB_Cell=cell(size(ImgStack_WS_Lb_Input,3),1);
Img_Gd_Lb_NoEdgeMapRGB_Cell=cell(size(ImgStack_WS_Lb_Input,3),1);


for k=1:size(ImgStack_WS_Lb_Input,3)

        % === Assign cells, bone / vessel / ImgMask DQ regions,
        % regionprops DQ regions to other colors, first by color index 
        % 0 values are edges of regions; do not erase

        % Watershed regions, Grid also identified as cells by "FindMaxima"
        Maxima_Cell_WSHullID=Maxima_Cell_WSHullID_Cell_Input{k,1};    
        Img_WS_Lb_IdxMap_Maxima=ImgStack_WS_Lb_Input(:,:,k);    
        Img_WS_Lb_IdxMap_Maxima(ismember(Img_WS_Lb_IdxMap_Maxima,Maxima_Cell_WSHullID(:,2)))=-2;
        
        % Grid identified as Watershed cells (with or without FindMaxima)
         Img_Gd_Lb_IdxMap_Cell=Img_Gd_Lb_Input;
         Img_Gd_Lb_IdxMap_Cell(ismember(Img_Gd_Lb_IdxMap_Cell,Gd_CellDQ_FlagID_Cell_Input{k,1}))=-3;
         
        % Watershed regions disqualified by region props
        Img_WS_Lb_IdxMap_PropDQ=ImgStack_WS_Lb_Input(:,:,k);   
        Img_WS_Lb_IdxMap_PropDQ(ismember(Img_WS_Lb_IdxMap_PropDQ,WS_Prop_FlagID_Cell_Input{k,1}))=-4;
         
        % Watershed regions, Grid disqualified by bone mask=bone area
        Img_WS_Lb_IdxMap_BoneDQ=ImgStack_WS_Lb_Input(:,:,k);   
        Img_WS_Lb_IdxMap_BoneDQ(ismember(Img_WS_Lb_IdxMap_BoneDQ,WS_BoneDQ_FlagID_Cell_Input{k,1}))=-6;
        Img_Gd_Lb_IdxMap_BoneDQ=Img_Gd_Lb_Input;
        Img_Gd_Lb_IdxMap_BoneDQ(ismember(Img_Gd_Lb_IdxMap_BoneDQ,Gd_BoneDQ_FlagID_Cell_Input{k,1}))=-6;
      
        % Watershed regions, Grid disqualified by vessel mask=vessel area 
        Img_WS_Lb_IdxMap_VesselDQ=ImgStack_WS_Lb_Input(:,:,k);   
        Img_WS_Lb_IdxMap_VesselDQ(ismember(Img_WS_Lb_IdxMap_VesselDQ,WS_VesselDQ_FlagID_Cell_Input{k,1}))=-7;  
        Img_Gd_Lb_IdxMap_VesselDQ=Img_Gd_Lb_Input;
        Img_Gd_Lb_IdxMap_VesselDQ(ismember(Img_Gd_Lb_IdxMap_VesselDQ,Gd_VesselDQ_FlagID_Cell_Input{k,1}))=-7;
        
        % Watershed regions, Grid disqualified by BCMask 
        Img_WS_Lb_IdxMap_BCMaskDQ=ImgStack_WS_Lb_Input(:,:,k);   
        Img_WS_Lb_IdxMap_BCMaskDQ(ismember(Img_WS_Lb_IdxMap_BCMaskDQ,WS_BCMaskDQ_FlagID_Cell_Input{k,1}))=-8;  
        Img_Gd_Lb_IdxMap_BCMaskDQ=Img_Gd_Lb_Input;
        Img_Gd_Lb_IdxMap_BCMaskDQ(ismember(Img_Gd_Lb_IdxMap_BCMaskDQ,Gd_BCMaskDQ_FlagID_Cell_Input{k,1}))=-8;
        
        % Watershed regions, Grid disqualified by ImgMask 
        Img_WS_Lb_IdxMap_ImgMaskDQ=ImgStack_WS_Lb_Input(:,:,k);   
        Img_WS_Lb_IdxMap_ImgMaskDQ(ismember(Img_WS_Lb_IdxMap_ImgMaskDQ,WS_ImgMaskDQ_FlagID_Cell_Input{k,1}))=-9;  
        Img_Gd_Lb_IdxMap_ImgMaskDQ=Img_Gd_Lb_Input;
        Img_Gd_Lb_IdxMap_ImgMaskDQ(ismember(Img_Gd_Lb_IdxMap_ImgMaskDQ,Gd_ImgMaskDQ_FlagID_Cell_Input{k,1}))=-9;
        
        % === Color Index, Watershed Labeled Map
        % Do not change order of color index assignment; preference implied
        % in coloring order
        Img_WS_Lb_IdxMap=ImgStack_WS_Lb_Input(:,:,k);
        
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap>0)=1;    % Use later; Watershed Cells *NOT* identified as cells by FindMaxima
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_Maxima==-2)=2; % Watershed Cells identified as cells by FindMaxima
        
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_PropDQ==-4)=4; % Watershed regions rejected by regionprops
        
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_BoneDQ==-6)=6; % Watershed regions rejected by bone mask
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_VesselDQ==-7)=7; % Watershed regions rejected by vessel mask
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_BCMaskDQ==-8)=8; % Watershed regions rejected by BCMask
        Img_WS_Lb_IdxMap(Img_WS_Lb_IdxMap_ImgMaskDQ==-9)=9; % Watershed regions rejected by ImgMask  
        
        % === Color Index, Grid Labeled Map
        % Do not change order of color index assignment; preference implied
        % in coloring order
        Img_Gd_Lb_IdxMap=single(logical(Img_Gd_Lb_Input)).*5; % Grid unassigned to cells / bone / vessels / ImgMask
        
        Img_Gd_Lb_IdxMap(Img_Gd_Lb_IdxMap_Cell==-3)=3; % Grid overlapping watershed Cells, both identified as cells by FindMaxima and not
        
        Img_Gd_Lb_IdxMap(Img_Gd_Lb_IdxMap_BoneDQ==-6)=6; % Grid overlapping watershed regions rejected by bone mask
        Img_Gd_Lb_IdxMap(Img_Gd_Lb_IdxMap_VesselDQ==-7)=7; % Grid overlapping watershed regions rejected by vessel mask
        Img_Gd_Lb_IdxMap(Img_Gd_Lb_IdxMap_BCMaskDQ==-8)=8; % Grid overlapping watershed regions rejected by BCMask
        Img_Gd_Lb_IdxMap(Img_Gd_Lb_IdxMap_ImgMaskDQ==-9)=9; % Grid overlapping watershed regions rejected by ImgMask
        
        % === Color Edges in Grid Labeled Map
        % Edge assigned to most prevalent neighboring feature
        % (cell/bone/vessel etc) w/ function "mode"
        EdgeRemv_SE = strel('square',3);    % Dilation for easier to read pts
        Img_Gd_Lb_NoEdgeIdxMap=Img_Gd_Lb_IdxMap;
        Img_Gd_Lb_NoEdgeIdxMap(Img_Gd_Lb_IdxMap==0)=NaN;
        Img_Gd_Lb_NoEdgeIdxMap=colfilt(Img_Gd_Lb_NoEdgeIdxMap,[3 3],'sliding',@mode);  % mode function ignores NaN
        Img_Gd_Lb_NoEdgeIdxMap(logical(Img_Gd_Lb_IdxMap))=-9999;  % Set Non-Edge pixel values to v. small
        Img_Gd_Lb_NoEdgeIdxMap=max(Img_Gd_Lb_IdxMap,Img_Gd_Lb_NoEdgeIdxMap); % Replace Non-Edge pixel values w/ Edge Map values
                    
        % === Assign actual colors by color index        
        % label2rgb output is uint8 class
        % Color assignment as follows:
        % 1) Watershed cells, not FindMaxima: Dull Yellow [0.8,0.8,0]
        % 2) Watershed cells, FindMaxima: Bright Yellow [1,1,0]
        % 3) Watershed cells, FindMaxima or not: Mid Yellow [0.9,0.9,0]
        % 4) Regionprops DQ regions: Grey [0.25,0.25,0.25]
        % 5) Not cells/bone/vessel/ImgMask (includes 4)): Dull Grey [0.125,0.125,0.125]
        % 6) Bone DQ regions: Blue [0,0,1]
        % 7) Vessel DQ regions: Red [1,0,0]
        % 8) BCMask DQ regions: Dark Green [1,0,1]
        % 9) ImgMask DQ regions: Dark Green [0,0.25,0]
        Lb_MapRGB_Color=...
            [0.80,0.80,0;...
            1,1,0;...
            0.90,0.90,0;...
            0.25,0.25,0.25;...
            0.125,0.125,0.125;...
            0,0,1;...
            1,0,0;...
            1,0,1;
            0,0.25,0];
        
        Img_WS_Lb_MapRGB=uint8(label2rgb(Img_WS_Lb_IdxMap,Lb_MapRGB_Color));    % label2rgb produced uint8 images
        Img_Gd_Lb_MapRGB=uint8(label2rgb(Img_Gd_Lb_IdxMap,Lb_MapRGB_Color));
        Img_Gd_Lb_NoEdgeMapRGB=uint8(label2rgb(Img_Gd_Lb_NoEdgeIdxMap,Lb_MapRGB_Color));
        
        % === Output
        ImgStack_WS_Lb_IdxMap(:,:,k)=Img_WS_Lb_IdxMap;
        ImgStack_Gd_Lb_IdxMap(:,:,k)=Img_Gd_Lb_IdxMap;
        ImgStack_Gd_Lb_NoEdgeIdxMap(:,:,k)=Img_Gd_Lb_NoEdgeIdxMap;
        
        Img_WS_Lb_MapRGB_Cell{k,1}=Img_WS_Lb_MapRGB;
        Img_Gd_Lb_MapRGB_Cell{k,1}=Img_Gd_Lb_MapRGB;
        Img_Gd_Lb_NoEdgeMapRGB_Cell{k,1}=Img_Gd_Lb_NoEdgeMapRGB;
    
end

clearvars ImgStack_WS_Lb;
clearvars Img_WS_Lb_IdxMap Img_Gd_Lb_IdxMap Img_Gd_Lb_NoEdgeIdxMap;
clearvars Img_WS_Lb_MapRGB Img_Gd_Lb_MapRGB Img_Gd_Lb_NoEdgeMapRGB;
clearvars Img_WS_Lb_IdxMap_Maxima;
clearvars Maxima_Cell_WSHullID;
clearvars Img_Gd_Lb_IdxMap_Cell;
clearvars Img_Gd_Lb_MapRGB_Cell;
clearvars Gd*FlagID_Cell;
clearvars Img_Gd_Lb;
clearvars *DQ*;
clearvars *_Input;


%% Colocalization check: Cell, bone, vessel, ImgMask as regions on Input Img (ImgStack)

fprintf('Colocalization check (ImgStack): Cells as regions ...\n');

Img_WS_Lb_MapRGB_Cell_Input=Img_WS_Lb_MapRGB_Cell;
Img_Gd_Lb_NoEdgeMapRGB_Cell_Input=Img_Gd_Lb_NoEdgeMapRGB_Cell;

WSMap_Img_Colocal_Cell=cell(NumImgSlices,1); % Output
GdNEMap_Img_Colocal_Cell=cell(NumImgSlices,1); % Output

% === Colocalization with ImgStack

for MapType=1:2 % Watershed region map, Grid map
        
        for k=1:NumImgSlices

                % Prepare colocalization layers
                if MapType==1   % Img_WS_Lb_IdxMapRGB
                    MapRGB_Layer=vertcat(zeros(Img_Height,Img_Width,3),Img_WS_Lb_MapRGB_Cell_Input{k,1});
                elseif MapType==2   % Img_Gd_Lb_NoEdgeMapRGB
                    MapRGB_Layer=vertcat(zeros(Img_Height,Img_Width,3),Img_Gd_Lb_NoEdgeMapRGB_Cell_Input{k,1});
                end
                Im_Layer=repmat(ImgStack(:,:,k).*255,2,1);
                Im_Layer=Im_Layer.*1/ImgStack_BCAdjMean;   % Adjust intensity for visualization
                Im_Layer(Im_Layer>255)=255;
                
                % Colocalize Map Layer and Img Layer
                FigHandle=figure('Toolbar','none','Menubar','none','Visible', 'off');
                imshow(uint8(Im_Layer),'InitialMagnification',80);
                hold on;
                MapHandle=imshow(uint8(MapRGB_Layer),'InitialMagnification',80);
                hold off;
                alpha=0.25;
                set(MapHandle, 'AlphaData', alpha);
                Fig=getframe(gcf);
                Map_Img_Colocal=Fig.cdata;
                close(FigHandle);
                
                if MapType==1   % Img_WS_Lb_IdxMapRGB
                    WSMap_Img_Colocal_Cell{k,1}=Map_Img_Colocal;
                elseif MapType==2   % Img_Gd_Lb_NoEdgeMapRGB
                    GdNEMap_Img_Colocal_Cell{k,1}=Map_Img_Colocal;
                end

        end
        
end

clearvars MapType;
clearvars MapRGB_Layer Im_Layer alpha Fig Map_Img_Colocal;
clearvars *_Input;


%% Save Images

fprintf('Saving Cell/Bone/Vessel - ImgStack Colocal Images...\n');

for k=1:NumImgSlices

    % Watershed Region Color Map of Cell, Bone, Vessel, ImgMask  
    % imwrite(Img_WS_Lb_MapRGB_Cell{k,1},strcat(SaveWSMapFilePath,'WSMap -',num2str(k,'%04.0f'),'.tif'));

    % Grid Color Map of Cell, Bone, Vessel, ImgMask 
    imwrite(Img_Gd_Lb_NoEdgeMapRGB_Cell{k,1},strcat(SaveGdMapFilePath,'GdMap -',num2str(k,'%04.0f'),'.tif'));
    
    % Colocalization: Watershed Region Color Map of Cell, Bone, Vessel, ImgMask  w/ ImgStack
    % imwrite(WSMap_Img_Colocal_Cell{k,1},strcat(SaveColocalWSMapImgFilePath,'ColocalWSMapImg -',num2str(k,'%04.0f'),'.tif'));
    
    % Colocalization: Grid Color Map of Cell, Bone, Vessel, ImgMask  w/ ImgStack
    % imwrite(GdNEMap_Img_Colocal_Cell{k,1},strcat(SaveColocalGdMapImgFilePath,'ColocalGdMapImg -',num2str(k,'%04.0f'),'.tif'));

end

clearvars WSMap_Img_Colocal_Cell GdNEMap_Img_Colocal_Cell;


%% Colocalization check: Cell, bone, vessel, ImgMask as regions on OrgImgStack
   
fprintf('Colocalization check (OrgImgStack; enhance contrast): Cells as regions ...\n');

Img_WS_Lb_MapRGB_Cell_Input=Img_WS_Lb_MapRGB_Cell;
Img_Gd_Lb_NoEdgeMapRGB_Cell_Input=Img_Gd_Lb_NoEdgeMapRGB_Cell;

WSMap_OrgImg_Colocal_Cell=cell(NumImgSlices,1); % Output
GdNEMap_OrgImg_Colocal_Cell=cell(NumImgSlices,1); % Output

% === Colocalization with ImgStack

for MapType=1:2 % Watershed region map, Grid map

    for k=1:NumImgSlices

            % Prepare colocalization layers
            if MapType==1   % Img_WS_Lb_IdxMapRGB
                MapRGB_Layer=vertcat(zeros(Img_Height,Img_Width,3),Img_WS_Lb_MapRGB_Cell_Input{k,1});
            elseif MapType==2   % Img_Gd_Lb_NoEdgeMapRGB
                MapRGB_Layer=vertcat(zeros(Img_Height,Img_Width,3),Img_Gd_Lb_NoEdgeMapRGB_Cell_Input{k,1});
            end
            Im_Layer=imadjust(OrgImgStack(:,:,k),stretchlim(OrgImgStack(:,:,k)));
            Im_Layer=repmat(Im_Layer.*255,2,1);
            Im_Layer=Im_Layer.*1;   % Adjust intensity for visualization
            Im_Layer(Im_Layer>255)=255; 

            % Colocalize Map Layer and Img Layer
            FigHandle=figure('Toolbar','none','Menubar','none','Visible', 'off');
            imshow(uint8(Im_Layer),'InitialMagnification',80);
            hold on;
            MapHandle=imshow(uint8(MapRGB_Layer),'InitialMagnification',80);
            hold off;
            alpha=0.25;
            set(MapHandle, 'AlphaData', alpha);
            Fig=getframe(gcf);
            Map_Img_Colocal=Fig.cdata;
            close(FigHandle);

            if MapType==1   % Img_WS_Lb_IdxMapRGB
                WSMap_OrgImg_Colocal_Cell{k,1}=Map_Img_Colocal;
            elseif MapType==2   % Img_Gd_Lb_NoEdgeMapRGB
                GdNEMap_OrgImg_Colocal_Cell{k,1}=Map_Img_Colocal;
            end

    end

end

clearvars MapType;
clearvars MapRGB_Layer Im_Layer alpha Fig Map_Img_Colocal;
clearvars *_Input;


%% Save Images

fprintf('Saving Cell/Bone/Vessel - OrgImgStack Colocal Images...\n');

for k=1:NumImgSlices

    % Colocalization: Watershed Region Color Map of Cell, Bone, Vessel, ImgMask  w/ OrgImgStack
    imwrite(WSMap_OrgImg_Colocal_Cell{k,1},strcat(SaveColocalWSMapOrgImgFilePath,'ColocalWSMapOrgImg -',num2str(k,'%04.0f'),'.tif'));

     % Colocalization: Grid Color Map of Cell, Bone, Vessel, ImgMask  w/ OrgImgStack
    % imwrite(GdNEMap_OrgImg_Colocal_Cell{k,1},strcat(SaveColocalGdMapOrgImgFilePath,'ColocalGdMapOrgImg -',num2str(k,'%04.0f'),'.tif'));

end

clearvars WSMap_OrgImg_Colocal_Cell GdNEMap_OrgImg_Colocal_Cell;
clearvars Img_WS_Lb_MapRGB_Cell Img_Gd_Lb_NoEdgeMapRGB_Cell;


%% Prepare Cell Boundary Map for colocalization check

fprintf('Cell boundary map for colocalization check ...\n');

% INPUT: Index Map of WS Cell location
ImgStack_WS_Lb_IdxMap_Input=ImgStack_WS_Lb_IdxMap;

ImgStack_WSNotMaximaCellBound=zeros(size(ImgStack_WS_Lb_IdxMap_Input));
ImgStack_WSMaximaCellBound=zeros(size(ImgStack_WS_Lb_IdxMap_Input));

for k=1:size(ImgStack_WS_Lb_IdxMap_Input,3)
% for k=1:1
    
    % Cells designated by watershed  but not by FindMaxima (after bone/vessel/ImgMask + regionprops DQ)
    Img_WSNotMaximaCellBound=zeros(Img_Height,Img_Width,'logical');
    
    WSNotMaximaCell=logical(ImgStack_WS_Lb_IdxMap_Input(:,:,k)==1);
    [WSNotMaximaCellBound,~]=bwboundaries(WSNotMaximaCell,'noholes');    
    WSNotMaximaCellBound_RowCol_List=cell2mat(WSNotMaximaCellBound);    % Combine boundary data of all cells
    if (~isempty(WSNotMaximaCellBound_RowCol_List))
        WSNotMaximaCellBound_LinIdx_List=sub2ind([Img_Height Img_Width],WSNotMaximaCellBound_RowCol_List(:,1),WSNotMaximaCellBound_RowCol_List(:,2)); % List of [Row,Col] of boundary
        Img_WSNotMaximaCellBound(WSNotMaximaCellBound_LinIdx_List)=0.5; 
    end
    
    % Cells designated by watershed AND FindMaxima (after bone/vessel/ImgMask + regionprops DQ)
    Img_WSMaximaCellBound=zeros(Img_Height,Img_Width,'logical');
    
    WSMaximaCell=logical(ImgStack_WS_Lb_IdxMap_Input(:,:,k)==2);
    [WSMaximaCellBound,~]=bwboundaries(WSMaximaCell,'noholes');    
    WSMaximaCellBound_RowCol_List=cell2mat(WSMaximaCellBound);    % Combine boundary data of all cells
    if (~isempty(WSMaximaCellBound_RowCol_List))
        WSMaximaCellBound_LinIdx_List=sub2ind([Img_Height Img_Width],WSMaximaCellBound_RowCol_List(:,1),WSMaximaCellBound_RowCol_List(:,2)); % List of [Row,Col] of boundary
        Img_WSMaximaCellBound(WSMaximaCellBound_LinIdx_List)=0.5;
    end
    
    % Output
    ImgStack_WSNotMaximaCellBound(:,:,k)=Img_WSNotMaximaCellBound;
    ImgStack_WSMaximaCellBound(:,:,k)=Img_WSMaximaCellBound;
        
end

clearvars WSNotMaximaCellBound WSNotMaximaCellBound_LinIdx_List WSNotMaximaCellBound_RowCol_List;
clearvars Img_WSNotMaximaCellBound Img_WSMaximaCellBound;
clearvars WSMaximaCell* WSNotMaximaCell*;
clearvars *_Input;


%% Colocalization check: Cell as boundaries on Input Img (ImgStack)

fprintf('Colocalization check (ImgStack): Cells as boundaries ...\n');

ImgStack_WSNotMaximaCellBound_Input=ImgStack_WSNotMaximaCellBound;
ImgStack_WSMaximaCellBound_Input=ImgStack_WSMaximaCellBound;

WSCellBound_Img_Colocal_Cell=cell(size(ImgStack_WSNotMaximaCellBound_Input,3),1); % Output

% === Colocalization with ImgStack

for k=1:size(ImgStack_WSNotMaximaCellBound_Input,3)
        
        WSCellBound_Img_Colocal=zeros(Img_Height*2,Img_Width,3,'single');
        
        % Prepare colocalization layers
        WSNoMaximaCellBound_Layer=vertcat(zeros(Img_Height,Img_Width,1),single(ImgStack_WSNotMaximaCellBound_Input(:,:,k)));
        WSMaximaCellBound_Layer=vertcat(zeros(Img_Height,Img_Width,1),single(ImgStack_WSMaximaCellBound_Input(:,:,k)));
        Im_Layer=repmat(ImgStack(:,:,k),2,1);
        Im_Layer=Im_Layer.*1/ImgStack_BCAdjMean;   % Adjust intensity for visualization
        Im_Layer(Im_Layer>1)=1; 

        % Colocalize Map Layer and Img Layer
        WSCellBound_Img_Colocal(:,:,1)=Im_Layer;   %Red
        WSCellBound_Img_Colocal(:,:,2)=WSMaximaCellBound_Layer;   %Green
        WSCellBound_Img_Colocal(:,:,3)=WSNoMaximaCellBound_Layer; %Blue
        
        % Output (type 'single')
        WSCellBound_Img_Colocal_Cell{k,1}=WSCellBound_Img_Colocal;

end

clearvars ImgStack;
clearvars ImgStack_BCAdjMean;
clearvars WSCellBound_Img_Colocal;
clearvars *_Layer;
clearvars *_Input;


%% Colocalization check: Cell as boundaries on Original Img Stack (OrgImgStack)

fprintf('Colocalization check (OrgImgStack; enhanced contrast): Cells as boundaries ...\n');

ImgStack_WSNotMaximaCellBound_Input=ImgStack_WSNotMaximaCellBound;
ImgStack_WSMaximaCellBound_Input=ImgStack_WSMaximaCellBound;

WSCellBound_OrgImg_Colocal_Cell=cell(size(ImgStack_WSNotMaximaCellBound_Input,3),1); % Output

% === Colocalization with OrgImgStack

for k=1:size(ImgStack_WSNotMaximaCellBound_Input,3)

        WSCellBound_OrgImg_Colocal=zeros(Img_Height*2,Img_Width,3,'single');

        % Prepare colocalization layers
        WSNoMaximaCellBound_Layer=vertcat(zeros(Img_Height,Img_Width,1),single(ImgStack_WSNotMaximaCellBound_Input(:,:,k)));
        WSMaximaCellBound_Layer=vertcat(zeros(Img_Height,Img_Width,1),single(ImgStack_WSMaximaCellBound_Input(:,:,k)));
        Im_Layer=imadjust(OrgImgStack(:,:,k),stretchlim(OrgImgStack(:,:,k)));
        Im_Layer=repmat(Im_Layer,2,1);
        Im_Layer=Im_Layer.*1;    % Adjust intensity for visualization
        Im_Layer(Im_Layer>1)=1; 

        % Colocalize Map Layer and Img Layer
        WSCellBound_OrgImg_Colocal(:,:,1)=Im_Layer;   % Red
        WSCellBound_OrgImg_Colocal(:,:,2)=WSMaximaCellBound_Layer;   % Green
        WSCellBound_OrgImg_Colocal(:,:,3)=WSNoMaximaCellBound_Layer; % Blue

        % Output (type 'single')
        WSCellBound_OrgImg_Colocal_Cell{k,1}=WSCellBound_OrgImg_Colocal;

end

clearvars OrgImgStack;
clearvars WSCellBound_OrgImg_Colocal;
clearvars *_Layer;
clearvars ImgStack_WSNotMaximaCellBound ImgStack_WSMaximaCellBound;
clearvars *_Input;


%% Save Images
% Colocalization: Cell Boundaries w/ ImgStack & OrgImgStack

fprintf('Saving Cell Boundary Colocalization Images...\n');

for k=1:NumImgSlices
        
    % imwrite(double(WSCellBound_Img_Colocal_Cell{k,1}),strcat(SaveColocalCellBoundImgFilePath,'ColocalCellBoundImg -',num2str(k,'%04.0f'),'.tif'));   
    imwrite(double(WSCellBound_OrgImg_Colocal_Cell{k,1}),strcat(SaveColocalCellBoundOrgImgFilePath,'ColocalCellBoundOrgImg -',num2str(k,'%04.0f'),'.tif'));
   
end

clearvars WSCellBound_Img_Colocal_Cell WSCellBound_OrgImg_Colocal_Cell;


%% Calculate Cell Count, Cell Area, Cell Density

fprintf('Calculating cell count, area and density ...\n');

Cell_Centroid_LinIdx_List_Cell_Input=WS_Cell_Centroid_LinIdx_List_Cell;  
WS_Cell_Area_List_Cell_Input=WS_Cell_Area_List_Cell;

Maxima_Cell_LinIdx_List_Cell_Input=Maxima_Cell_LinIdx_List_Cell;
Maxima_Cell_Area_List_Cell_Input=Maxima_Cell_Area_List_Cell;

ImgStack_WS_Lb_IdxMap_Input=ImgStack_WS_Lb_IdxMap;
ImgStack_Gd_Lb_NoEdgeIdxMap_Input=ImgStack_Gd_Lb_NoEdgeIdxMap;

% Output
CellCtAreaDensity_metric_Mtx=zeros(size(ImgStack_WS_Lb_IdxMap_Input,3),23);
CellCtAreaDensity_pxUnit_Mtx=zeros(size(ImgStack_WS_Lb_IdxMap_Input,3),23);

for k=1:NumImgSlices
    
        WS_Cell_Centroid_LinIdx_List=Cell_Centroid_LinIdx_List_Cell_Input{k,1};
        WS_Cell_Area_List=WS_Cell_Area_List_Cell_Input{k,1};
        Maxima_Cell_LinIdx_List=Maxima_Cell_LinIdx_List_Cell_Input{k,1};
        Maxima_Cell_Area_List=Maxima_Cell_Area_List_Cell_Input{k,1};
        Img_WS_Lb_IdxMap=ImgStack_WS_Lb_IdxMap_Input(:,:,k);
        Img_Gd_Lb_NoEdgeIdxMap=ImgStack_Gd_Lb_NoEdgeIdxMap_Input(:,:,k);

        % === Cell Count (exclude bone/vessel/ImgMask, regionprops)

        % Watershed 
        if ~isempty(WS_Cell_Centroid_LinIdx_List)
            CellCt_WS=size(WS_Cell_Centroid_LinIdx_List,1);  
        else
            CellCt_WS=0;
        end

        % FIndMaxima
        if Maxima_Cell_LinIdx_List_Cell_Input{k,1}(1,1)==-9999
            CellCt_Maxima=0;
        else
            CellCt_Maxima=size(Maxima_Cell_LinIdx_List,1);   %Cell Count by FindMaxima (after regionprop, env reject)
        end

        % === Area, Diameter
  
        % Bone (Gd map)
        Area_GdBone_px2=length(find(Img_Gd_Lb_NoEdgeIdxMap==6));
        Area_GdBone_um2=Area_GdBone_px2*(XY_PxLength^2);

        % Vessel (Gd map)
        Area_GdVessel_px2=length(find(Img_Gd_Lb_NoEdgeIdxMap==7));
        Area_GdVessel_um2=Area_GdVessel_px2*(XY_PxLength^2);

        % BCMask (Gd map)
        Area_GdBCMask_px2=length(find(Img_Gd_Lb_NoEdgeIdxMap==8));
        Area_GdBCMask_um2=Area_GdBCMask_px2*(XY_PxLength^2);
        
        % ImgMask (Gd map)
        Area_GdImgMask_px2=length(find(Img_Gd_Lb_NoEdgeIdxMap==9));
        Area_GdImgMask_um2=Area_GdImgMask_px2*(XY_PxLength^2);

        % Cavity (includes area occupied by vessels)
        Area_GdCavity_px2=Img_Height*Img_Width-Area_GdBone_px2-Area_GdBCMask_px2-Area_GdImgMask_px2;
        Area_GdCavity_um2=Area_GdCavity_px2*(XY_PxLength^2);

        % Marrow Space (excludes area occupied by vessels)
        Area_GdMarrow_px2=Img_Height*Img_Width-Area_GdBone_px2-Area_GdVessel_px2-Area_GdBCMask_px2-Area_GdImgMask_px2;
        Area_GdMarrow_um2=Area_GdMarrow_px2*(XY_PxLength^2);        
        
        % All Cells (WS map; NOTE: includes leakage areas)
        if ~isempty(WS_Cell_Area_List)
            Area_CellWS_px2=sum(WS_Cell_Area_List,1);
            Area_Mean_CellWS_px2=mean(WS_Cell_Area_List,1);
            Area_Stdev_CellWS_px2=std(WS_Cell_Area_List,0,1);
        else
            Area_CellWS_px2=0;
            Area_Mean_CellWS_px2=0;
            Area_Stdev_CellWS_px2=0;
        end
         
        Area_CellWS_um2=Area_CellWS_px2*(XY_PxLength^2);
        Area_Mean_CellWS_um2=Area_Mean_CellWS_px2*(XY_PxLength^2);
        Area_Stdev_CellWS_um2=Area_Stdev_CellWS_px2*(XY_PxLength^2);
        
        if ~isempty(WS_Cell_Area_List)
            WS_Cell_Diameter_List=sqrt(WS_Cell_Area_List/pi)*2;
            Diameter_Mean_CellWS_px=mean(WS_Cell_Diameter_List,1);
            Diameter_Stdev_CellWS_px=std(WS_Cell_Diameter_List,0,1);
        else
            WS_Cell_Diameter_List=0;
            Diameter_Mean_CellWS_px=0;
            Diameter_Stdev_CellWS_px=0;
        end
        
        Diameter_Mean_CellWS_um=Diameter_Mean_CellWS_px*XY_PxLength;
        Diameter_Stdev_CellWS_um=Diameter_Stdev_CellWS_px*XY_PxLength;
              
        % Maxima Cells Only (WS map; NOTE: includes leakage areas)
        if ~isempty(Maxima_Cell_Area_List)
            if Maxima_Cell_Area_List(1)~=-9999 
                Area_CellMaxima_px2=sum(Maxima_Cell_Area_List,1);   
                Area_Mean_CellMaxima_px2=mean(Maxima_Cell_Area_List,1);   
                Area_Stdev_CellMaxima_px2=std(Maxima_Cell_Area_List,0,1);   
            else
                Area_CellMaxima_px2=0;   
                Area_Mean_CellMaxima_px2=0; 
                Area_Stdev_CellMaxima_px2=0;
            end
        else
                Area_CellMaxima_px2=0;   
                Area_Mean_CellMaxima_px2=0; 
                Area_Stdev_CellMaxima_px2=0;
        end
               
        Area_CellMaxima_um2=Area_CellMaxima_px2*(XY_PxLength^2);  
        Area_Mean_CellMaxima_um2=Area_Mean_CellMaxima_px2*(XY_PxLength^2);  
        Area_Stdev_CellMaxima_um2=Area_Stdev_CellMaxima_px2*(XY_PxLength^2);  
        
        if ~isempty(Maxima_Cell_Area_List)
            if Maxima_Cell_Area_List(1)~=-9999
                Maxima_Cell_Diameter_List=sqrt(Maxima_Cell_Area_List/pi)*2;
                Diameter_Mean_CellMaxima_px=mean(Maxima_Cell_Diameter_List,1);
                Diameter_Stdev_CellMaxima_px=std(Maxima_Cell_Diameter_List,0,1);
            else
                Maxima_Cell_Diameter_List=0;
                Diameter_Mean_CellMaxima_px=0;
                Diameter_Stdev_CellMaxima_px=0;
            end
        else
            Maxima_Cell_Diameter_List=0;
            Diameter_Mean_CellMaxima_px=0;
            Diameter_Stdev_CellMaxima_px=0;
        end
        
        Diameter_Mean_CellMaxima_um=Diameter_Mean_CellMaxima_px*XY_PxLength;
        Diameter_Stdev_CellMaxima_um=Diameter_Stdev_CellMaxima_px*XY_PxLength;
                
        % === Density 
        % Cell Ct/ Marrow Area (excludes vessel area)
        if Area_GdMarrow_px2 > 0
            Density_CellCtWS_Marrow_px2=CellCt_WS/Area_GdMarrow_px2;
            Density_CellCtWS_Marrow_mm2=CellCt_WS/Area_GdMarrow_um2*1E6;
            
            Density_CellCtMaxima_Marrow_px2=CellCt_Maxima/Area_GdMarrow_px2;
            Density_CellCtMaxima_Marrow_mm2=CellCt_Maxima/Area_GdMarrow_um2*1E6;
        else
            Density_CellCtWS_Marrow_px2=-9999;
            Density_CellCtWS_Marrow_mm2=-9999;
            
            Density_CellCtMaxima_Marrow_px2=-9999;
            Density_CellCtMaxima_Marrow_mm2=-9999;
        end

         % Cell Ct/ Cavity Area (includes vessel area)
         if Area_GdCavity_px2 > 0
            Density_CellCtWS_Cavity_px2=CellCt_WS/Area_GdCavity_px2;
            Density_CellCtWS_Cavity_mm2=CellCt_WS/Area_GdCavity_um2*1E6;
             
            Density_CellCtMaxima_Cavity_px2=CellCt_Maxima/Area_GdCavity_px2;
            Density_CellCtMaxima_Cavity_mm2=CellCt_Maxima/Area_GdCavity_um2*1E6;            
         else
            Density_CellCtWS_Cavity_px2=-9999;
            Density_CellCtWS_Cavity_mm2=-9999;
             
            Density_CellCtMaxima_Cavity_px2=-9999;
            Density_CellCtMaxima_Cavity_mm2=-9999;
         end

        % === Output
        % Load Cell Ct Area Density Data into Data Matrix;
        % Separate into 2 data files: 1 for metric unit measurements, 
        % one for pixel unit measurement
        
        % Data matrix contain info for cell ct, area, cell density for all slices

        % Col 01: SliceNum
        % Col 02: Cell Ct, by Watershed (WS) centroid (after bone/vessel/ImgMask+regionprop reject)
        % Col 03: Cell Ct, by ImgJ FindMaxima (Maxima) (after bone/vessel/ImgMask+regionprop reject)

        % Col 04: Area of Marrow (ImgSize-ImgMaskArea-BoneArea-VesselArea)
        % Col 05: WS Cell density in Marrow (Col 02/Col 04)
        % Col 06: Maxima Cell density in Marrow (Col 03/Col 04)

        % Col 07: Area of Cavity (ImgSize-ImgMaskArea-BoneArea)
        % Col 08: WS Cell density in Cavity (Col 02/Col 07)
        % Col 09: Maxima Cell density in Marrow (Col 03/Col 07)

        % Col 10: Area of Bone
        % Col 11: Area of Vessel
        % Col 12: Area of BCMask
        % Col 13: Area of ImgMask

        % Col 14: Total Area of (all) Watershed Cells in slice
        % Col 15: Mean Area of (each) Watershed Cell
        % Col 16: Standard Deviation Area of (each) Watershed Cell
        % Col 17: Mean Diameter of (each) Watershed Cell
        % Col 18: Standard Deviation Diameter of (each) Watershed Cell
        
        % Col 19: Total Area of (all) Maxima Cells in slice
        % Col 20: Mean Area of (each) Maxima Cell
        % Col 21: Standard Deviation Area of (each) Maxima Cell
        % Col 22: Mean Diameter of (each) Maxima Cell
        % Col 23: Standard Deviation Diameter of (each) Maxima Cell
        
        % = Metric
        CellCtAreaDensity_metric_Mtx(k,:)=...
            [k,...          
            CellCt_WS,...
            CellCt_Maxima,...
            Area_GdMarrow_um2,...
            Density_CellCtWS_Marrow_mm2,...
            Density_CellCtMaxima_Marrow_mm2,...
            Area_GdCavity_um2,...
            Density_CellCtWS_Cavity_mm2,...
            Density_CellCtMaxima_Cavity_mm2,...
            Area_GdBone_um2,...
            Area_GdVessel_um2,...
            Area_GdBCMask_um2,...
            Area_GdImgMask_um2,...
            Area_CellWS_um2,...
            Area_Mean_CellWS_um2,...
            Area_Stdev_CellWS_um2,...
            Diameter_Mean_CellWS_um,...
            Diameter_Stdev_CellWS_um,...
            Area_CellMaxima_um2,...
            Area_Mean_CellMaxima_um2,...
            Area_Stdev_CellMaxima_um2,...
            Diameter_Mean_CellMaxima_um,...
            Diameter_Stdev_CellMaxima_um];
        
        % = Pixel unit
        CellCtAreaDensity_pxUnit_Mtx(k,:)=...
            [k,...            
            CellCt_WS,...
            CellCt_Maxima,...
            Area_GdMarrow_px2,...
            Density_CellCtWS_Marrow_px2,...
            Density_CellCtMaxima_Marrow_px2,...
            Area_GdCavity_px2,...
            Density_CellCtWS_Cavity_px2,...
            Density_CellCtMaxima_Cavity_px2,...
            Area_GdBone_px2,...
            Area_GdVessel_px2,...
            Area_GdBCMask_px2,...
            Area_GdImgMask_px2,...
            Area_CellWS_px2,...
            Area_Mean_CellWS_px2,...
            Area_Stdev_CellWS_px2,...
            Diameter_Mean_CellWS_px,...
            Diameter_Stdev_CellWS_px,...
            Area_CellMaxima_px2,...
            Area_Mean_CellMaxima_px2,...
            Area_Stdev_CellMaxima_px2,...
            Diameter_Mean_CellMaxima_px,...
            Diameter_Stdev_CellMaxima_px];
        
end

% Display approximate Cell Density results on screen
WSCellDensity_Display=CellCtAreaDensity_metric_Mtx(:,5);
WSCellDensity_Display=WSCellDensity_Display(WSCellDensity_Display>eps);
MaximaCellDensity_Display=CellCtAreaDensity_metric_Mtx(:,6);
MaximaCellDensity_Display=MaximaCellDensity_Display(MaximaCellDensity_Display>eps);
fprintf('\n====================\n');
fprintf(strcat('WS Cell density in marrow space (mm2)~',num2str(mean(WSCellDensity_Display(:,1))),'\n'));
% fprintf(strcat('Maxima Cell density in marrow space (mm2)~',num2str(mean(MaximaCellDensity_Display(:,1))),'\n'));
fprintf('=====================\n\n');

clearvars WS_Cell_Centroid_LinIdx_List WS_Cell_Area_List WS_Cell_Diameter_List;
clearvars Maxima_Cell_LinIdx_List Maxima_Cell_Area_List Maxima_Cell_Diameter_List;
clearvars ImgStack_WS_Lb_IdxMap ImgStack_Gd_Lb_IdxMap ImgStack_Gd_Lb_NoEdgeIdxMap;
clearvars Img_WS_Lb_IdxMap Img_Gd_Lb_IdxMap Img_Gd_Lb_NoEdgeIdxMap;
clearvars Area_* Diameter_* Density_* CellCt_*;

clearvars *_Input;


%% Histogram, Cell Diameter

fprintf('Cell Diameter Histogram...\n');

WS_Cell_Area_List_Cell_Input=WS_Cell_Area_List_Cell;
Maxima_Cell_Area_List_Cell_Input=Maxima_Cell_Area_List_Cell;

% Calculate cell diameter 
WS_Cell_Area_List_AllSlice_um2=cell2mat(WS_Cell_Area_List_Cell_Input)*(XY_PxLength^2);
WS_Cell_Diameter_List_AllSlice_um=sqrt(WS_Cell_Area_List_AllSlice_um2/pi)*2;

Maxima_Cell_Area_List_AllSlice_um2=cell2mat(Maxima_Cell_Area_List_Cell_Input)*(XY_PxLength^2);
Maxima_Cell_Area_List_AllSlice_um2(Maxima_Cell_Area_List_AllSlice_um2<0)=0;    % Remove -9999 entries
Maxima_Cell_Diameter_List_AllSlice_um=sqrt(Maxima_Cell_Area_List_AllSlice_um2/pi)*2;

% Display approximate Cell Diameter results on screen
fprintf('\n=====================\n');
fprintf(strcat('WS Cell mean diameter (um)=',num2str(mean(WS_Cell_Diameter_List_AllSlice_um,1)),'\n'));
fprintf(strcat('WS Cell stdev diameter (um)=',num2str(std(WS_Cell_Diameter_List_AllSlice_um,0,1)),'\n'));
% fprintf(strcat('Maxima Cell mean diameter (um)=',num2str(mean(Maxima_Cell_Diameter_List_AllSlice_um,1)),'\n'));
% fprintf(strcat('Maxima Cell stdev diameter (um)=',num2str(std(Maxima_Cell_Diameter_List_AllSlice_um,0,1)),'\n'));
fprintf('=====================\n\n');

% Plot and save Cell Diameter histograms
figHandle01=figure;
histogram(WS_Cell_Diameter_List_AllSlice_um,100);
ylabel('Histogram Ct'); 
xlabel('WS Cell Diameter, AllStack (um)');
title({strcat('WS Cell mean diameter (um)=',num2str(mean(WS_Cell_Diameter_List_AllSlice_um,1))),...
    strcat('WS Cell stdev diameter (um)=',num2str(std(WS_Cell_Diameter_List_AllSlice_um,0,1))),...
    strcat('WS Cell density in marrow space (mm2)~',num2str(mean(WSCellDensity_Display(:,1))))});
print(figHandle01,'-dtiffn','-r0',strcat(SaveFilePath,'WSCellSizeHistogram_',timestamp,'.tif'));

% figHandle02=figure;
% histogram(Maxima_Cell_Diameter_List_AllSlice_um,100);
% ylabel('Histogram Ct'); 
% xlabel('Maxima Cell Diameter, AllStack (um)');
% title({strcat('Maxima Cell mean diameter (um)=',num2str(mean(Maxima_Cell_Diameter_List_AllSlice_um,1))),...
%     strcat('Maxima Cell stdev diameter (um)=',num2str(std(Maxima_Cell_Diameter_List_AllSlice_um,0,1))),...
%     strcat('Maxima Cell density in marrow space (mm2)~',num2str(mean(MaximaCellDensity_Display(:,1))))});
% print(figHandle02,'-dtiffn','-r0',strcat(SaveFilePath,'MaximaCellSizeHistogram_',timestamp,'.tif'));

% Close all figures
close all;

clearvars *_Display;
clearvars *AllSlice*;
clearvars *_Input;


%% Save cell count, area, density data

fprintf('Saving cell count, area and density data...\n');

CellCtAreaDensity_metric_Mtx_Input=CellCtAreaDensity_metric_Mtx;
CellCtAreaDensity_pxUnit_Mtx_Input=CellCtAreaDensity_pxUnit_Mtx;

for UNITTYPE=1:2 
    
        if UNITTYPE==1    % Metric units
            Data_Mtx=CellCtAreaDensity_metric_Mtx_Input;
            Data_FileNameString=strcat(SaveFilePath,'CellCtAreaDensityMetric AllSliceData.txt');
            Data_FileHeaderRow={'SliceNum';...
                'CellCt,WatershedWS';'CellCt,FindMaximaFM';...
                'MarrowArea(um2)';'MarrowCellDensity,WS(mm-2)';'MarrowCellDensity,FM(mm-2)';...
                'CavityArea(um2)';'CavityCellDensity,WS(mm-2)';'CavityCellDensity,FM(mm-2)';...
                'BoneArea(um2)';'VesselArea(um2)';'BCMaskArea(um2)';'ImgMaskArea(um2)';...
                'TotalCellArea,WS(um2)';'MeanCellArea,WS(um2)';'StdevCellArea,WS(um2)';'MeanCellDiameter,WS(um)';'StdevCellDiameter,WS(um)';...
                'TotalCellArea,FM(um2)';'MeanCellArea,FM(um2)';'StdevCellArea,FM(um2)';'MeanCellDiameter,FM(um)';'StdevCellDiameter,FM(um)'};
        elseif UNITTYPE==2    % Pixel units
            Data_Mtx=CellCtAreaDensity_pxUnit_Mtx;
            Data_FileNameString=strcat(SaveFilePath,'CellCtAreaDensityPxUnit AllSliceData.txt');
            Data_FileHeaderRow={'SliceNum';...
                'CellCt,WatershedWS';'CellCt,FindMaximaFM';...
                'MarrowArea(px2)';'MarrowCellDensity,WS(px-2)';'MarrowCellDensity,FM(px-2)';...
                'CavityArea(px2)';'CavityCellDensity,WS(px-2)';'CavityCellDensity,FM(px-2)';...
                'BoneArea(px2)';'VesselArea(px2)';'BCMaskArea(px2)';'ImgMaskArea(px2)';...
                'TotalCellArea,WS(px2)';'MeanCellArea,WS(px2)';'StdevCellArea,WS(px2)';'MeanCellDiameter,WS(px)';'StdevCellDiameter,WS(px)';...
                'TotalCellArea,FM(px2)';'MeanCellArea,FM(px2)';'StdevCellArea,FM(px2)';'MeanCellDiameter,FM(px)';'StdevCellDiameter,FM(px)'};        
        end

        % = Save
        Data_FileID=fopen(Data_FileNameString,'w');

        % Write Header rows
        for i=1:numel(Data_FileHeaderRow)
            fprintf(Data_FileID,'%s\t',Data_FileHeaderRow{i});
        end
        fprintf(Data_FileID,'\n');

        % Write data
        dlmwrite(Data_FileNameString,Data_Mtx,'delimiter','\t','-append');

        fclose(Data_FileID);

end

clearvars Data_*;
clearvars *_Input;


%% Save Watershed Cell Position Coordinate + Area
% Save X-coordinate, Y-coordinate and Linear Index of Watershed Cell
% Centroids

fprintf('Saving cell position coordinates and cell area...\n');

WS_Cell_Centroid_XY_List_Cell_Input=WS_Cell_Centroid_XY_List_Cell;      
WS_Cell_Centroid_LinIdx_List_Cell_Input=WS_Cell_Centroid_LinIdx_List_Cell;    
WS_Cell_Area_List_Cell_Input=WS_Cell_Area_List_Cell;

Maxima_Cell_XY_List_Cell_Input=Maxima_Cell_XY_List_Cell;      
Maxima_Cell_LinIdx_List_Cell_Input=Maxima_Cell_LinIdx_List_Cell;  
Maxima_Cell_Area_List_Cell_Input=Maxima_Cell_Area_List_Cell;

% for METHOD=1:2
for METHOD=1:1
    
        if METHOD==1    % Watershed Centroid 
            CellPos_XY_List_Cell_Input=WS_Cell_Centroid_XY_List_Cell_Input;
            CellPos_LinIdx_List_Cell_Input=WS_Cell_Centroid_LinIdx_List_Cell_Input; 
            CellArea_Cell_Input=WS_Cell_Area_List_Cell_Input; 
            CellPosArea_FileNameString=strcat(SaveFilePath,'WatershedCellPosArea AllSliceData.txt');
            CellPosArea_FileHeaderRow={'SliceNum';'WSCellPos,X';'WSCellPos,Y';'WSCellPos,LinIdx';'WSCellArea(px2)'};
        elseif METHOD==2    % Maxima
            CellPos_XY_List_Cell_Input= Maxima_Cell_XY_List_Cell_Input;
            CellPos_LinIdx_List_Cell_Input=Maxima_Cell_LinIdx_List_Cell_Input;
            CellArea_Cell_Input=Maxima_Cell_Area_List_Cell_Input; 
            CellPosArea_FileNameString=strcat(SaveFilePath,'MaximaCellPosArea AllSliceData.txt');
            CellPosArea_FileHeaderRow={'SliceNum';'MaximaCellPos,X';'MaximaCellPos,Y';'MaximaCellPos,LinIdx';'WSCellArea(px2)'};
        end
        
        % === Create Data Matrix
        CellPosArea_Mtx=zeros(size(cell2mat(CellPos_XY_List_Cell_Input),1),5);
        
        % k Index for Matrix loading
        CellPosArea_Mtx_RowIdx=0;   

        for k=1:size(CellPos_XY_List_Cell_Input,1)

            SliceNum=repmat(k,size(CellPos_XY_List_Cell_Input{k,1},1),1);    
            CellPos_XY=CellPos_XY_List_Cell_Input{k,1};
            CellPos_LinIdx=CellPos_LinIdx_List_Cell_Input{k,1};
            CellArea=CellArea_Cell_Input{k,1};
            
            if ~isempty(CellPos_XY)
                CellPosArea=horzcat(SliceNum,CellPos_XY,CellPos_LinIdx,CellArea);
            else
                CellPosArea=zeros(0,5); % Empty row
            end
            
            CellPosArea_Mtx(CellPosArea_Mtx_RowIdx+1:CellPosArea_Mtx_RowIdx+size(CellPosArea,1),:)=CellPosArea;    
            CellPosArea_Mtx_RowIdx=CellPosArea_Mtx_RowIdx+size(CellPosArea,1);

        end
        
        % Delete rows with -9999 values
        CellPosArea_Mtx(CellPosArea_Mtx(:,4)==-9999,:)=[];

        % === Save
        CellPosArea_FileID=fopen(CellPosArea_FileNameString,'w');

        % Write header row
        for i=1:numel(CellPosArea_FileHeaderRow)
            fprintf(CellPosArea_FileID,'%s\t',CellPosArea_FileHeaderRow{i});
        end
        fprintf(CellPosArea_FileID,'\n');

        % Write data
        dlmwrite(CellPosArea_FileNameString,CellPosArea_Mtx,'delimiter','\t','-append','precision', 16);

        fclose(CellPosArea_FileID);

end

clearvars WS_Cell_*_Cell;
clearvars Maxima_Cell_*_Cell;
clearvars Maxima_NoWSHull_LinIdx_List_Cell Maxima_Cell_WSHullID_Cell;
clearvars CellPos*;
clearvars *_Input;

clearvars Save*Path -except SaveFilePath;
clearvars j k i regionct gridct ItCt SliceNum METHOD UNITTYPE scrsz timestamp;
clearvars Img;
clearvars CellArea;
clearvars *Handle*;


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);


%% Save workspace

% clearvars IJM MIJ;
%  
% Workspace_FileNameString=strcat(SaveFilePath,'CellCtWorkspace',timestamp_str,'.mat');
% save(Workspace_FileNameString);

