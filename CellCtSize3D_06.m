clc; fclose('all'); clearvars;
close all hidden;

addpath(genpath('INPUT/')) % Add INPUT folder and subfolders to path
addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path


%% =====DESCRIPTION=====

% Estimate 3D cell density and size distribution based on 2D cell centroid, size and position info

% == Usage:
% User specifies variables in "USER INPUT" section.

% ==Output folders:
% (Each folder includes .MAT file and its corresponding 8-bit image stack)
% "CellCentroid3D_OrgImgStack_Colocal_ZSAR"? Co-localized 3D centroids and original vessel channel image stack

% ==Output files:
% "PostCutCellDiameter3DHistogram.tif": 3D cell size distribution histogram


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA
% Email address: jwwu@@mgh.harvard.edu
% Last revision: June-2021


%% USER INPUT

% === INPUT: ImgStack (GdMap from CellCtSize2D.m)
% InputImg Tif vs MAT Query:
% TIF: img from microscope (ex. 1st pass)
% MAT: Single channel single precision matrix
InputImg_TifvsMAT_Query=1;

% Input image directory (do NOT include / at end)
CellCtSize2DOut_Folder_struct=dir(fullfile('OUTPUT/','CellCtSize2D*'));
BatchImgInputFolder=strcat('OUTPUT/',CellCtSize2DOut_Folder_struct(end).name,'/GdMap');

% Input images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255)
ImgStackData_FilenameString='';


% === INPUT: CellPosArea AllSlice (txt file from CellCtSize2D.m)
CellPosArea_FilenameString=strcat('OUTPUT/',CellCtSize2DOut_Folder_struct(end).name,'/WatershedCellPosArea AllSliceData.txt');


% === INPUT: Original image colocalization
% 3D Cell Positions (Centroid) is overlaid on original images with brightness and contrast
% adjustments only
OrgImgColocal_Query='y';

% Input original image directory (do NOT include / aImg_Widtht end)
InputOrgImg_TifvsMAT_Query=2;

% Input original image directory (do NOT include / at end)
BatchOrgImgInputFolder='';

% Input original images: if .MAT, code reads image intensity matrix (in greyscale [0,1]*255)
OrgImgStackData_FilenameString='INPUT/DEMO_ImgStack_G.mat';


% === INPUT: Output image directory (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/CellCtSize3D_Out';


% === INPUT: Pixel length in XY, in Z (um)
XY_PxLength=0.3405;
Z_PxLength=1.0;


% === INPUT: Centroid Position Match Query
% 'y' if match 3D Centroid to one of 2D's
% 'n; if do not match (use weighted centroid)
CentroidPosMatch_Query='y';


% === INPUT: Mininum Cell Diameter, um
% TargetMinCellDiameter_um is an estimate
% Actual value depends on ThreshCellRadius_Ratio and HistEdge
% If supplying CellCentroid3DSizePosAllCells_FilenameString,
% Check MinCellDiameter_um_LoLimit
TargetMinCellDiameter_um=4.2;


% === INPUT: Reduce processing time by limiting ZSlice to process
ManualSliceRange_Query='n';
ManualStartSlice=1;
ManualEndSlice=107;


% === INPUT: Limiting Z Range for Cells
% Does not change bone / vessel / marrow cavity volumes.
ManualCellSliceRange_Query='y';
ManualCellStartSlice=1; % Slice included
ManualCellEndSlice=30;


%% Optimized inputs: modify with care

% === INPUT: ThreshCellRadius_px_Ratio
ThreshCellRadius_Ratio=0.5; % Do not change

% === INPUT: Stepping for Radius, in pixels
CellRadius_px_Step_Query=0.2;

% === INPUT: Dimension of cell position gradient map
% Determine if 2D or 3D Gaussian is used to re-generate each "cell"
Dim_CellPosGradMap_Query=3; % 2 or 3

% === INPUT: for imgaussfilt(3)
% = Diameter/Sigma Ratio
DiameterSigmaRatio_Query=5;

% = Multiplier for imgaussfilt(3)
% Dimension for Multiplier for CellPosGradMap, based on Radius
% 0 = Dim 0; multiplier = 1
% 1 = Dim 1; multiplier = Diameter
% 2 = Dim 2 (2D analysis) or 3 (3D analysis); multiplier = area (2D) or
% volume (3D)
GradMapMultDim_Max_Query=1; % Do not change
GradMapMultDim_Sum_Query=2; % Do not change

% === INPUT: Option to load CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx
CellCentroid3DSizePosAllCells_Query='n';

% Name of CellCentroid3DSizePosAllCells File
CellCentroid3DSizePosAllCells_FilenameString='';

% === INPUT: Gdmap Color
% These colors should match those specified in CellCtSize2D.m
BoneColor=[0,0,1]; % [R,G,B]
VesselColor=[1,0,0]; % [R,G,B]
BCMaskColor=[1,0,1]; % [R,G,B]
ImgMaskColor=[0,0.25,0]; % [R,G,B]

% === INPUT: Color of OrgImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=2;


%% Load GdMap to get Img_Height, Img_Width and NumImgSlices info only

% Load image stack (GdMap)
[~,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad_f(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,1);


%% Read CellPosArea information
% Calculate 3D Linear Idx for Position And Approximate Diameter

% 1st column: Slice number
% 2nd column: Cell centroid position, Row number
% 3rd column: Cell centroid position, Column number
% 4rd column: Cell centroid position, LinIdx (by 2D plane)
% 5th column: Cell area, in pixel^2

fprintf('Reading Cell Position and Area Data...\n');

CellPosArea_FileID=fopen(CellPosArea_FilenameString);
CellPosArea=dlmread(CellPosArea_FilenameString,'\t',1,0);

if ManualCellSliceRange_Query=='y'
    CellPosArea(find(CellPosArea(:,1)<ManualCellStartSlice),:)=[];
    CellPosArea(find(CellPosArea(:,1)>ManualCellEndSlice),:)=[];
end


%% Option to limit to 250 random cells for testing

% % % fprintf('Limiting to 250 cells (for debugging)...\n');
% % %
% % % Randvector=round(rand(250,1)*size(CellPosArea,1));
% % % Randvector(Randvector==0)=[];
% % % CellPosArea=CellPosArea(Randvector,:);


%% Prepare for Saving

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/');
mkdir(SaveFilePath);


%% Determine # Z slices for processing
% Lowest and Highest Z for cells

fprintf('Preparing for Z shift correction..\n');

CellPosArea_Input=CellPosArea;

if ManualSliceRange_Query=='y'
    StartImgSliceZS=ManualStartSlice;
    EndImgSliceZS=ManualEndSlice;
    if StartImgSliceZS>min(CellPosArea_Input(:,1))
        fprintf('\nWARNING: StartImgSlice_Cells below 1st ZSlice containing 2D cell centroids...\n\n');
    end
    NumImgSlicesZS=ManualEndSlice-ManualStartSlice+1; % Pre Aspect Ratio Fix
    
else
    StartImgSliceZS=1;
    EndImgSliceZS=NumImgSlices;
    NumImgSlicesZS=NumImgSlices; % Pre Aspect Ratio Fix
end

clearvars *_Input;


%% Z-Shift: Shift ZSlice of 2D Cell Centroid ; set 1:1:1 X:Y:Z aspect ratio
% Voxel length in X,Y,Z set to XY_PxLength, adjust number of Z slices accordingly

fprintf('Performing: ...\n');
fprintf('1) Z-Shift (ZS): Shift starting ZSlice to ManualStartSlice, if selected.\n');
fprintf('2) Aspect Ratio Fix (AR): Converting cell position Data to 1:1:1 X:Y:Z aspect ratio.\n');
fprintf('Voxel size set to XY_PxLength, adjust number of Z slices accordingly...\n');

CellPosArea_Input=CellPosArea;

% Cut off unecessary rows before Z-shift
CellPosArea_Input(CellPosArea_Input(:,1)<StartImgSliceZS,:)=[];
CellPosArea_Input(CellPosArea_Input(:,1)>=EndImgSliceZS,:)=[];

% Z-Shift
CellPosAreaZS=CellPosArea_Input; %ZS= S Shift
CellPosAreaZS(:,1)=CellPosAreaZS(:,1)-StartImgSliceZS+1;

% Aspect Ratio Fix
CellPosAreaZSAR=CellPosAreaZS;
CellPosAreaZSAR(:,1)=round(CellPosAreaZSAR(:,1)*Z_PxLength/XY_PxLength);
NumImgSlicesZSAR=round(NumImgSlicesZS*Z_PxLength/XY_PxLength);
Z_PxLength_ZSAR=Z_PxLength*NumImgSlicesZS/NumImgSlicesZSAR;

clearvars CellPosArea CellPosAreaZS;
clearvars *_Input;


%% Create 2D cell diameter (um) histogram
% This should be identical to the histogram generated in Step 10

fprintf('Creating 2D cell diameter (um) histogram...\n');

CellPosArea_Input=CellPosAreaZSAR;

% = Plot 2D Cell Radius Data
% Calculate approximate diameter by assuming area is circle
% Used for both 2D and 3D
CellRadiusXY_px=sqrt(CellPosArea_Input(:,5)/pi);    % Pre-Pad

% Define Cell Radius Histogram Edge, will also be used later in script.
% All cells in the same histogram group will be processed together
% Group CellRadius by histogram, in steps of 0.2px.
CellRadius_px_Step_Query=0.2;
CellRadiusMin=max(CellRadius_px_Step_Query,(floor(min(CellRadiusXY_px(:,1)*round(1/CellRadius_px_Step_Query))))/round(1/CellRadius_px_Step_Query));
CellRadiusMax=(ceil(max(CellRadiusXY_px(:,1)*round(1/CellRadius_px_Step_Query))))/round(1/CellRadius_px_Step_Query);
CellRadius_HistEdge=(CellRadiusMin:CellRadius_px_Step_Query:CellRadiusMax);    % To be Used later
CellDiameter2D_HistEdge_um=CellRadius_HistEdge*XY_PxLength*2;

% FigHandle00=figure('Position',[10 300 600 600]);
% histogram(CellRadiusXY_px*XY_PxLength*2,CellDiameter2D_HistEdge_um,'FaceColor',[1 1 0]);
% ylabel('Histogram: 2D Cell Diameter (um)');
% xlabel('2D Cell Diameter (um)');
% title({strcat('2D Cell Diameter (um)'),...
%     strcat('Total Cell Count (2D)=',num2str(size(CellPosArea_Input,1)))});
% print(FigHandle00,'-dtiffn','-r0',strcat(SaveFilePath,'CellDiameter2DHistogram.tif'));
% close(FigHandle00);

clearvars CellRadiusXY_px;
clearvars *_Input;


%% Prepare CellCentroid3D Information

if CellCentroid3DSizePosAllCells_Query=='y'
    
    fprintf('Loading 3D cell count, centroid information...\n');
    
    CellCentroid3DSizePosAllCells=load(CellCentroid3DSizePosAllCells_FilenameString);
    
    % Check CellCentroid3DSizePosAllCells_FilenameString information
    % matches other input
    if ~isequal(CellPosAreaZSAR,CellCentroid3DSizePosAllCells.CellPosAreaZSAR)
        error('Error. CellPosAreaZSAR and CellPosAreaZSAR in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(CellRadius_HistEdge,CellCentroid3DSizePosAllCells.CellRadius_HistEdge)
        error('Error. CellRadius_HistEdge and CellRadius_HistEdge in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(CellRadius_px_Step_Query,CellCentroid3DSizePosAllCells.CellRadius_px_Step_Query)
        error('Error. CellRadius_px_Step_Query and CellRadius_px_Step_Query in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(DiameterSigmaRatio_Query,CellCentroid3DSizePosAllCells.DiameterSigmaRatio_Query)
        error('Error. DiameterSigmaRatio_Query and DiameterSigmaRatio_Query in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(ThreshCellRadius_Ratio,CellCentroid3DSizePosAllCells.ThreshCellRadius_Ratio)
        error('Error. ThreshCellRadius_Ratio and ThreshCellRadius_Ratio in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(GradMapMultDim_Max_Query,CellCentroid3DSizePosAllCells.GradMapMultDim_Max_Query)
        error('Error. GradMapMultDim_Max_Query and GradMapMultDim_Max_Query in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    if ~isequal(GradMapMultDim_Sum_Query,CellCentroid3DSizePosAllCells.GradMapMultDim_Sum_Query)
        error('Error. GradMapMultDim_Sum_Query and GradMapMultDim_Sum_Query in "CellCentroid3DSizePosAllCells_FilenameString" does not match.')
    end
    
    CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx=CellCentroid3DSizePosAllCells.CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx;
    MinCellDiameter_um_LoLimit=CellCentroid3DSizePosAllCells.MinCellDiameter_um_LoLimit;
    CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(find(CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(:,1)<MinCellDiameter_um_LoLimit),:)=[];
    
    
elseif CellCentroid3DSizePosAllCells_Query=='n'
    
    %% Create Gradient Stack using Cell Position and Area
    % Step 1: Prepare data matrix based on 2D cell radius grouping
    
    fprintf('Creating 3D cell gradient stack using 2D cell position and area data.\n');
    fprintf('Step 1: Prepare data matrix based on 2D cell radius grouping...\n');
    
    CellPosArea_Input=CellPosAreaZSAR;
    CellRadiusHistGrp_HistEdge_Input=CellRadius_HistEdge;
    
    % = Padding
    CellPosArea_PrePosPad=CellPosArea_Input;
    CellPosArea_PosPad=CellPosArea_Input;
    CellPosArea_PosPadStep=-round((Z_PxLength/XY_PxLength-1)/2):1:round((Z_PxLength/XY_PxLength-1)/2);  % Z_PxLength, NOT Z_PxLength_ZSAR
    for PadCt=1:size(CellPosArea_PosPadStep,2)
        if abs(CellPosArea_PosPadStep(PadCt))>eps   % Ignore 0
            CellPosArea_PosPadTemp=CellPosArea_PrePosPad;
            CellPosArea_PosPadTemp(:,1)=CellPosArea_PosPadTemp(:,1)+CellPosArea_PosPadStep(1,PadCt);
            CellPosArea_PosPad=vertcat(CellPosArea_Input,CellPosArea_PosPadTemp);
        end
    end
    
    % Re-Calculate approximate diameter by assuming area is circle
    % Used for both 2D and 3D
    CellRadiusXY_px_PosPad=sqrt(CellPosArea_PosPad(:,5)/pi);
    
    % Combine CellPos (in 3DLinIdx) and Radius info into 1 matrix for easy processing
    % Calculate 3D LinIdx of cell positions, referenced to StartImgSliceZS
    % (i.e. set StartImgSlice_Cells as Slice 1)
    % Used for both 2D and 3D
    CellPos_3DLinIdx_PosPad=CellPosArea_PosPad(:,1)*Img_Height*Img_Width+CellPosArea_PosPad(:,4);
    CellPos_Radius_Mtx=horzcat(CellPos_3DLinIdx_PosPad,CellRadiusXY_px_PosPad);
    
    % Group CellRadius by histogram, in steps of 0.2px.
    % All cells in the same histogram group will be processed together
    [CellRadiusHistGrp_HistCt,~,CellRadiusHistGrpIdx]=...
        histcounts(CellPos_Radius_Mtx(:,2),CellRadiusHistGrp_HistEdge_Input);
    NumCellRadiusHistGrp=numel(CellRadiusHistGrp_HistEdge_Input);
    
    % Split CellPos_Radius_Mtx into different matrices in cells according to
    % their CellRadius histogram group; allows parfor
    CellPos_Radius_RadiusHistGrp_Cell=cell(NumCellRadiusHistGrp,1);
    for HistGrpCt=1:NumCellRadiusHistGrp   % "for" faster than "parfor"
        CellPos_Radius_RadiusHistGrp_Cell{HistGrpCt,1}=CellPos_Radius_Mtx(CellRadiusHistGrpIdx==HistGrpCt,:);
    end
    
    clearvars CellPosArea_PrePosPad CellPosArea_PosPad CellPosArea_PosPadTemp CellPosArea_PosPadStep;
    clearvars PadCt;
    clearvars CellPos_3DLinIdx_PosPad;
    clearvars CellRadiusXY_px_PosPad CellPos_3DLinIdx_PosPad;
    clearvars CellRadiusHistGrpIdx;
    clearvars HistGrpCt;
    clearvars *_Input;
    
    
    %% Create Gradient Stack using Cell Position and Area
    % Step 2: Calculate 3D Gaussian peak intensity based on 2D cell radius grouping
    
    fprintf('Creating 3D cell gradient stack using 2D cell position and area data.\n');
    fprintf('Step 2: Calculate 3D Gaussian peak intensity based on 2D cell radius grouping...\n');
    
    % Calculate by applying imgaussfilt(3) to Temp stack
    % Use: deviate from theoretical values depending on CellRadius
    % as imgaussfilt(3) only allows positive, odd filter size
    Gauss_PDF_PeakIntensity_RadiusHistGrp=zeros(NumCellRadiusHistGrp,1);
    
    for  HistGrpCt=1:NumCellRadiusHistGrp
        CellRadiusXY=CellRadius_HistEdge(1,HistGrpCt);
        ImgStack_Temp_Size=DiameterSigmaRatio_Query*ceil(2*(CellRadiusXY/DiameterSigmaRatio_Query))+1;
        imgaussFiltSize=DiameterSigmaRatio_Query*ceil(2*(CellRadiusXY/DiameterSigmaRatio_Query))+1;
        if mod(imgaussFiltSize,2)>0
        else
            imgaussFiltSize=imgaussFiltSize+1;
        end
        if Dim_CellPosGradMap_Query==2
            ImgStack_Temp=zeros(ImgStack_Temp_Size,ImgStack_Temp_Size);
            ImgStack_Temp(round(ImgStack_Temp_Size/2),round(ImgStack_Temp_Size/2))=1;
            ImgStack_Temp=imgaussfilt(ImgStack_Temp,CellRadiusXY/DiameterSigmaRatio_Query,'FilterSize',imgaussFiltSize);
        elseif Dim_CellPosGradMap_Query==3
            ImgStack_Temp=zeros(ImgStack_Temp_Size,ImgStack_Temp_Size,ImgStack_Temp_Size);
            ImgStack_Temp(round(ImgStack_Temp_Size/2),round(ImgStack_Temp_Size/2),round(ImgStack_Temp_Size/2))=1;
            ImgStack_Temp=imgaussfilt3(ImgStack_Temp,CellRadiusXY/DiameterSigmaRatio_Query,'FilterSize',imgaussFiltSize);
        end
        Gauss_PDF_PeakIntensity_RadiusHistGrp(HistGrpCt,1)=max(ImgStack_Temp(:));
    end
    
    Gauss_PDF_PeakIntensity_RadiusHistGrp_Cell=...
        mat2cell(Gauss_PDF_PeakIntensity_RadiusHistGrp,ones(NumCellRadiusHistGrp,1));
    
    clearvars CellRadiusXY Gauss_coor_step;
    clearvars ImgStack_Temp ImgStack_Temp_Size;
    clearvars Gauss_sigma Gauss_sigma_covar Gauss_mu;
    clearvars MeshX MeshY MeshZ Gauss_coor;
    clearvars Gauss_PDF_3D Gauss_PDF_PeakIntensity_RadiusHistGrp;
    clearvars HistGrpCt;
    clearvars *_Input;
    
    
    %% Create Gradient Stack using Cell Position and Area
    % Step 3: Prepare image stack for 3D Gaussian Filtering with 2D cell data,
    % based on 2D cell radius grouping
    % Each 2D cell position is marked by point with intensity = area (2D) or volume (3D),
    % divided by 3D Gaussian peak intensity
    
    fprintf('Creating 3D cell gradient stack using 2D cell position and area data.\n');
    fprintf('Step 3: 3D Gaussian Filtering of 2D cell data, based on 2D cell radius grouping...\n');
    
    % To save memory, skip input line:
    % ImgStack_RadiusHistGrp_Cell_Input=ImgStack_RadiusHistGrp_Cell;
    imgaussFiltSize_Input=imgaussFiltSize;
    
    % GradMapMult indicates multiplier for each imgaussfilt(3) outcome stack
    % Multiplier choices based on dimension:
    % Unity (dim 0), Diameter (dim 1), Area (dim 2), Voume (dim 3)
    % Cell Column 1 for Max, Cell Column 2 for Sum
    GradMapMult_RadiusHistGrp_Cell=cell(NumCellRadiusHistGrp,2);
    
    for HistGrpCt=1:NumCellRadiusHistGrp
        
        % For CellPosGradMapMaxStack
        if GradMapMultDim_Max_Query==0
            GradMapMult_RadiusHistGrp_Cell{HistGrpCt,1}=1;
        elseif GradMapMultDim_Max_Query==1
            GradMapMult_RadiusHistGrp_Cell{HistGrpCt,1}=CellRadius_HistEdge(1,HistGrpCt)*2;
        elseif GradMapMultDim_Max_Query==2
            if Dim_CellPosGradMap_Query==2
                GradMapMult_RadiusHistGrp_Cell{HistGrpCt,1}=CellRadius_HistEdge(1,HistGrpCt).^2*pi;
            elseif Dim_CellPosGradMap_Query==3
                GradMapMult_RadiusHistGrp_Cell{HistGrpCt,1}=CellRadius_HistEdge(1,HistGrpCt).^3*4/3*pi;
            end
        end
        
        % For CellPosGradMapSumStack
        if GradMapMultDim_Sum_Query==0
            GradMapMult_RadiusHistGrp_Cell{HistGrpCt,2}=1;
        elseif GradMapMultDim_Sum_Query==1
            GradMapMult_RadiusHistGrp_Cell{HistGrpCt,2}=CellRadius_HistEdge(1,HistGrpCt)*2;
        elseif GradMapMultDim_Sum_Query==2
            if Dim_CellPosGradMap_Query==2
                GradMapMult_RadiusHistGrp_Cell{HistGrpCt,2}=CellRadius_HistEdge(1,HistGrpCt).^2*pi;
            elseif Dim_CellPosGradMap_Query==3
                GradMapMult_RadiusHistGrp_Cell{HistGrpCt,2}=CellRadius_HistEdge(1,HistGrpCt).^3*4/3*pi;
            end
        end
        
    end
    
    % ImgStack_RadiusHistGrp_Cell: one stack for each 2D cell radius grouping.
    % Each 2D cell position is marked by point with intensity = 1.
    ImgStack_RadiusHistGrp_Cell=cell(NumCellRadiusHistGrp,1);
    
    % CellPosGradMapMax_RadiusHistGrp_Cell: post 3D Gaussian filtering of
    % ImgStack_RadiusHistGrp_Cell
    CellPosGradMap_RadiusHistGrp_Cell=cell(NumCellRadiusHistGrp,1);
    
    % CellPosGradMapStack: assemble information from CellPosGradMap_RadiusHistGrp_Cell
    CellPosGradMapMaxStack=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
    CellPosGradMapSumStack=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
    
    tic
    parfor HistGrpCt=1:NumCellRadiusHistGrp  % "parfor" faster
        
        ImgStack_RadiusHistGrp_Cell{HistGrpCt,1}=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
        CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
        
        if ~isempty(CellPos_Radius_RadiusHistGrp_Cell{HistGrpCt,1})
            
            ImgStack_RadiusHistGrp_Cell{HistGrpCt,1}(CellPos_Radius_RadiusHistGrp_Cell{HistGrpCt,1}(:,1))=1;
            
            if Dim_CellPosGradMap_Query==2
                CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}=...
                    imgaussfilt(ImgStack_RadiusHistGrp_Cell{HistGrpCt,1},(CellRadius_HistEdge(1,HistGrpCt)/DiameterSigmaRatio_Query),'FilterSize',imgaussFiltSize_Input);
            elseif Dim_CellPosGradMap_Query==3
                CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}=...
                    imgaussfilt3(ImgStack_RadiusHistGrp_Cell{HistGrpCt,1},(CellRadius_HistEdge(1,HistGrpCt)/DiameterSigmaRatio_Query),'FilterSize',imgaussFiltSize_Input);
            end
            CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}=...
                CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}./Gauss_PDF_PeakIntensity_RadiusHistGrp_Cell{HistGrpCt,1};
            CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}(CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}>1)=1;
            
            CellPosGradMapMaxStack=...
                max(CellPosGradMapMaxStack,CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}.*GradMapMult_RadiusHistGrp_Cell{HistGrpCt,1}); % maximum
            CellPosGradMapSumStack=...
                CellPosGradMapSumStack+CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}.*GradMapMult_RadiusHistGrp_Cell{HistGrpCt,2}; % sum
            
            ImgStack_RadiusHistGrp_Cell{HistGrpCt,1}=[]; % save memory
            CellPosGradMap_RadiusHistGrp_Cell{HistGrpCt,1}=[]; % save memory
            
        end
        
    end
    toc
    
    clearvars ImgStack_RadiusHistGrp_Cell;
    clearvars CellPos_Radius_RadiusHistGrp_Cell;
    clearvars CellPosGradMap_RadiusHistGrp_Cell;
    clearvars GradMapMult_RadiusHistGrp_Cell;
    clearvars HistGrpCt;
    
    
    %% Display, Save Cell Position Gradient Map
    
    fprintf('Saving Cell Position Gradient Map...\n');
    
    if Dim_CellPosGradMap_Query==2
        SaveCellPosGradMapMaxStackFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/Gauss2D_Max_ZSAR',num2str(StartImgSliceZS),'-',num2str(EndImgSliceZS),'/');
        SaveCellPosGradMapSumStackFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/Gauss2D_Sum_ZSAR',num2str(StartImgSliceZS),'-',num2str(EndImgSliceZS),'/');
    elseif Dim_CellPosGradMap_Query==3
        SaveCellPosGradMapMaxStackFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/Gauss3D_Max_ZSAR',num2str(StartImgSliceZS),'-',num2str(EndImgSliceZS),'/');
        SaveCellPosGradMapSumStackFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/Gauss3D_Sum_ZSAR',num2str(StartImgSliceZS),'-',num2str(EndImgSliceZS),'/');
    end
    
    %         mkdir(SaveCellPosGradMapMaxStackFilePath);
    %         mkdir(SaveCellPosGradMapSumStackFilePath);
    
    % Adjust intensity to [0 1] for image saving only
    CellPosGradMapMaxStack_ImgSave=CellPosGradMapMaxStack;
    CellPosGradMapMaxStack_ImgSave=single((CellPosGradMapMaxStack_ImgSave-min(CellPosGradMapMaxStack_ImgSave(:)))*1/(max(CellPosGradMapMaxStack_ImgSave(:))-min(CellPosGradMapMaxStack_ImgSave(:))));
    
    CellPosGradMapSumStack_ImgSave=CellPosGradMapSumStack;
    CellPosGradMapSumStack_ImgSave=single((CellPosGradMapSumStack_ImgSave-min(CellPosGradMapSumStack_ImgSave(:)))*1/(max(CellPosGradMapSumStack_ImgSave(:))-min(CellPosGradMapSumStack_ImgSave(:))));
    
    %         figure;
    %         for k=1:size(CellPosGradMapMaxStack_ImgSave,3)
    %             imshow(CellPosGradMapMaxStack_ImgSave(:,:,k));
    %             imwrite(double(CellPosGradMapMaxStack_ImgSave(:,:,k)),...
    %                 strcat(SaveCellPosGradMapMaxStackFilePath,'CellPosGradMapMaxStackBCAdj -',num2str(k,'%04.0f'),'.tif'));
    %         end
    %
    %         for k=1:size(CellPosGradMapSumStack_ImgSave,3)
    %             imshow(CellPosGradMapSumStack_ImgSave(:,:,k));
    %             imwrite(double(CellPosGradMapSumStack_ImgSave(:,:,k)),...
    %                 strcat(SaveCellPosGradMapSumStackFilePath,'CellPosGradMapSumStackBCAdj -',num2str(k,'%04.0f'),'.tif'));
    %         end
    
    clearvars CellPosGradMapMaxStack_ImgSave CellPosGradMapSumStack_ImgSave;
    clearvars SaveCellPosGradMapMaxStackFilePath SaveCellPosGradMapSumStackFilePath;
    
    
    %% 3D cell count, locate centroid
    
    fprintf('Performing 3D cell count, locate centroid...\n');
    
    CellPosGradMapMaxStack_Input=CellPosGradMapMaxStack;
    CellPosGradMapSumStack_Input=CellPosGradMapSumStack;
    
    % Determine threshold level for normalized CellPosGradMapStack
    % Use XY_PxLength as scaling factor as sigma for Gaussians were
    % determined based on CellRadiusXY
    ThreshLevelScaleFactor=max(CellPosGradMapMaxStack_Input(:));
    
    CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx=[-9999,-9999]; % Dummy row
    
    % Create HistEdge for for-loop. At least 0.5px spacing.
    % Reverse order (largest radius first)
    CellCentroid3D_Radius_HistEdge=CellRadius_HistEdge(end:-ceil(0.5/CellRadius_px_Step_Query):1);
    MinCellDiameter_um_LoLimit=(CellCentroid3D_Radius_HistEdge(end-1)*2)/ThreshCellRadius_Ratio*XY_PxLength;
    
    for j=1:(size(CellCentroid3D_Radius_HistEdge,2)-1)
        
        % ThreshCellRadius_px=CellRadius_HistEdge(j)-CellRadius_px_Step_Query*0.5;
        ThreshCellRadius_px=CellCentroid3D_Radius_HistEdge(j)*ThreshCellRadius_Ratio;
        
        if ThreshCellRadius_px<CellCentroid3D_Radius_HistEdge(end)
            break;
        end
        
        if GradMapMultDim_Max_Query==0
            ThreshLevel=multiththresh(CellPosGradMapMaxStack_Input./ThreshLevelScaleFactor,1);
        elseif GradMapMultDim_Max_Query==1  % Diameter, in px
            ThreshLevel=ThreshCellRadius_px*2;
            ThreshLevel=ThreshLevel/ThreshLevelScaleFactor;
        elseif GradMapMultDim_Max_Query==2
            if Dim_CellPosGradMap_Query==2 % Area, in px
                ThreshLevel=(ThreshCellRadius_px).^2*pi;
                ThreshLevel=ThreshLevel/ThreshLevelScaleFactor;
            elseif Dim_CellPosGradMap_Query==3 % Volume, in px
                ThreshLevel=(ThreshCellRadius_px).^3*pi*4/3;
                ThreshLevel=ThreshLevel/ThreshLevelScaleFactor;
            end
        end
        
        % = Identify  & count cells in 3D, using 26-neighbour connectivity
        % Binarize CellPosGradMapMAXStack to get volumes to keep
        % Get 3D weighted centroid from CellPosGradMapSUMStack
        Cell3D_BW=imbinarize(CellPosGradMapMaxStack_Input/ThreshLevelScaleFactor,ThreshLevel);
        Cell3D_CC=bwconncomp(Cell3D_BW,26); % 3D connected components
        Cell3D_regionprops3=...
            regionprops3(Cell3D_CC,CellPosGradMapSumStack_Input/max(CellPosGradMapSumStack_Input(:)),'WeightedCentroid','VoxelIdxList');
        
        % = Determine largest radius covered by each 3D cell
        ImgStack_CellRadius_px=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
        
        if CentroidPosMatch_Query=='y'
            ImgStack_CellRadius_px((CellPosAreaZSAR(:,1)-1)*Img_Height*Img_Width+CellPosAreaZSAR(:,4))=sqrt(CellPosAreaZSAR(:,5)/pi);
        elseif CentroidPosMatch_Query=='n'
            ImgStack_CellRadius_px(CellPos_Radius_Mtx(:,1))=CellPos_Radius_Mtx(:,2); % Original
        end
        
        CellDiameter3D_px_Cell=cell(size(Cell3D_regionprops3,1),1);
        CellCentroid3D_VoxelIdx_Cell=cell(size(Cell3D_regionprops3,1),1);
        for CellCt=1:size(Cell3D_regionprops3,1) % "for" faster than parfor
                   [CellRadius3D_px_Temp,CellRadius3D_px_Temp_Idx]=max(ImgStack_CellRadius_px(Cell3D_regionprops3.VoxelIdxList{CellCt,1}));
            CellDiameter3D_px_Cell{CellCt,1}=CellRadius3D_px_Temp*2;
            CellCentroid3D_VoxelIdx_Cell{CellCt,1}=Cell3D_regionprops3.VoxelIdxList{CellCt,1}(CellRadius3D_px_Temp_Idx);
        end
        CellDiameter3D_um=cell2mat(CellDiameter3D_px_Cell)*(XY_PxLength*2+Z_PxLength_ZSAR)/3;
        CellCentroid3D_VoxelIdx=cell2mat(CellCentroid3D_VoxelIdx_Cell);
        
        % = Elimate 3D cells with Radius=0;
        fprintf(strcat('...',num2str(numel(find(CellDiameter3D_um<eps))),'/',num2str(size(Cell3D_regionprops3,1)),' cells to be removed due to lack of 3D radius information...\n'));
        Cell3D_regionprops3(CellDiameter3D_um<eps,:)=[];
        CellCentroid3D_VoxelIdx(CellDiameter3D_um<eps,:)=[];
        CellDiameter3D_um(CellDiameter3D_um<eps,:)=[];
        
        % = Elimate 3D cells withCell 3D_MaxDiameter_um not within
        % CellCentroid3D_Radius_Hist bin
        fprintf(strcat('...',num2str(numel(union(find(CellDiameter3D_um>=(CellCentroid3D_Radius_HistEdge(j)*2)*XY_PxLength),find(CellDiameter3D_um<(CellCentroid3D_Radius_HistEdge(j+1)*2)*XY_PxLength)))),'/',num2str(size(Cell3D_regionprops3,1)),' cells to be removed for diameter out of bounds...\n'));
        Cell3D_regionprops3(CellDiameter3D_um>=(CellCentroid3D_Radius_HistEdge(j)*2)*XY_PxLength,:)=[];
        CellCentroid3D_VoxelIdx(CellDiameter3D_um>=(CellCentroid3D_Radius_HistEdge(j)*2)*XY_PxLength,:)=[];
        CellDiameter3D_um(CellDiameter3D_um>=(CellCentroid3D_Radius_HistEdge(j)*2)*XY_PxLength,:)=[];
        Cell3D_regionprops3(CellDiameter3D_um<(CellCentroid3D_Radius_HistEdge(j+1)*2)*XY_PxLength,:)=[];
        CellCentroid3D_VoxelIdx(CellDiameter3D_um<(CellCentroid3D_Radius_HistEdge(j+1)*2)*XY_PxLength,:)=[];
        CellDiameter3D_um(CellDiameter3D_um<(CellCentroid3D_Radius_HistEdge(j+1)*2)*XY_PxLength,:)=[];
        
        if ~isempty(CellDiameter3D_um)
            
            % = Determine Centroid Values in Linear Index, Plot
            % Centroid value in [x,y,z]=[Col#,Row#,Slice#]
            
            if CentroidPosMatch_Query=='y'
                CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx=vertcat(CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx,horzcat(CellDiameter3D_um,CellCentroid3D_VoxelIdx));
            elseif CentroidPosMatch_Query=='n'
                CellCentroid3D_LinIdx=sub2ind(size(CellPosGradMapMaxStack_Input),...
                    round(Cell3D_regionprops3.WeightedCentroid(:,2)),round(Cell3D_regionprops3.WeightedCentroid(:,1)),round(Cell3D_regionprops3.WeightedCentroid(:,3)));
                CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx=vertcat(CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx,horzcat(CellDiameter3D_um,CellCentroid3D_LinIdx));
            end
            
        end
        
    end
    
    % Delete 1st dummy row, remove cells < MinCellDiameter_um_LoLimit
    CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(1,:)=[];
    CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(find(CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(:,1)<MinCellDiameter_um_LoLimit),:)=[];
    
    % Save CellCentroid3D_Diameter_um_LinIdx_Mtx before truncating based on MinCellDiameter_um
    %         save(strcat(SaveFilePath,'CellCentroid3DSizePosAllCells.mat'),...
    %             'CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx',...
    %             'CellRadius_HistEdge','CellPosAreaZSAR',...
    %             'CellRadius_px_Step_Query','DiameterSigmaRatio_Query',...
    %             'GradMapMultDim_Max_Query','GradMapMultDim_Sum_Query',...
    %             'ThreshCellRadius_Ratio','MinCellDiameter_um_LoLimit',...
    %             'ManualSliceRange_Query','ManualStartSlice','ManualEndSlice',...
    %             'ManualCellSliceRange_Query','ManualCellStartSlice','ManualCellEndSlice',...
    %             'Img_Height','Img_Width','NumImgSlicesZSAR');
   
    
end

clearvars CellDiameter3D_um CellCentroid3D_LinIdx;
clearvars CellCentroid3D_Radius_HistEdge;
clearvars ThreshLevel ThreshLevelScaleFactor;
clearvars ThreshCellRadius_px;
clearvars CellRadiusHistGrp_HistCt;
clearvars CellRadiusHistGrp_HistCt_MaxIdx;
clearvars CellPos_Radius_Mtx;
clearvars CellPosGradMapMaxStack CellPosGradMapSumStack;
clearvars Cell3D_BW Cell3D_CC Cell3D_LabelMtx;
clearvars ImgStack_CellRadius_px;
clearvars Cell3D_regionprops3;
clearvars CellDiameter3D_px_Cell;
clearvars CellCt;
clearvars CellPosAreaZSAR;
clearvars *_Input;


%% Create 3D cell diameter (um) histogram
% This should be identical to the histogram generated in Step 10

fprintf('Creating 3D cell diameter (um) histogram...\n');

CellDiameter3D_um_Input=CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx(:,1);

% For 3D cell diameter histogram, round off to nearest and larger 0.1
% based on the value of CellRadius_px_Step.
% (Rounding off here allows easier comparison in future)
CellDiameter3D_um_Step=ceil(CellRadius_px_Step_Query/((XY_PxLength*2+Z_PxLength_ZSAR)/3)*10)/10;
CellDiameter3D_um_HistEdge_min=floor(min(CellDiameter3D_um_Input)/CellDiameter3D_um_Step)*CellDiameter3D_um_Step;
CellDiameter3D_um_HistEdge_max=ceil(max(CellDiameter3D_um_Input)/CellDiameter3D_um_Step)*CellDiameter3D_um_Step;
CellDiameter3D_um_HistEdge=...
    CellDiameter3D_um_HistEdge_min:CellDiameter3D_um_Step:CellDiameter3D_um_HistEdge_max;
CellDiameter3D_um_HistCtBinLoc=0.5*(CellDiameter3D_um_HistEdge(1:end-1)+CellDiameter3D_um_HistEdge(2:end));
[CellDiameter3D_um_HistCt,~]=histcounts(CellDiameter3D_um_Input,CellDiameter3D_um_HistEdge);

% Set MinCellDiameter to coincide with smaller HistEdge value closest to TargetMinCellDiameter_um
% (So that the 1st bar is the correct visual representation of histogram bin)
% Also must not be < MinCellDiameter_um_LoLimit
if TargetMinCellDiameter_um<MinCellDiameter_um_LoLimit
    TargetMinCellDiameter_um=MinCellDiameter_um_LoLimit;
end
Diff=TargetMinCellDiameter_um-CellDiameter3D_um_HistEdge;
Diff(Diff<0)=9999; % TargetMinCellDiameter_um must > CellDiameter3D_um_HistEdge
[~,MinIdx]=min(Diff);
MinCellDiameter_um=CellDiameter3D_um_HistEdge(MinIdx);
while MinCellDiameter_um<MinCellDiameter_um_LoLimit
    MinIdx=MinIdx+1;
    MinCellDiameter_um=CellDiameter3D_um_HistEdge(MinIdx);
end
fprintf('\n... Based on TargetMinCellDiameter_um and method settings...\n');
fprintf(strcat('... MinCellDiameter_um=',num2str(MinCellDiameter_um),'.\n\n'));

% FigHandle01a=figure('Position',[30 300 600 600]);
% % histogram(CellDiameter3D_um_Input,CellDiameter3D_um_HistEdge,'FaceColor',[1 0.75 0]);
% bar(CellDiameter3D_um_HistCtBinLoc,CellDiameter3D_um_HistCt,'FaceColor',[1 0.75 0],'BarWidth',1);
% ylabel('Histogram: 3D Cell Diameter (um)');
% xlabel('3D Cell Diameter (um)');
% xticks(CellDiameter3D_um_HistCtBinLoc);
% xtickangle(45);
% ax=gca;
% ax.FontSize=18;
% title({strcat('3D Cell Diameter (um)'),...
%     strcat('Total Cell Count (3D)=',num2str(size(CellDiameter3D_um_Input,1))),...
%     strcat('Radius/Sigma for Gaussian=',num2str(DiameterSigmaRatio_Query)),...
%     strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
% print(FigHandle01a,'-dtiffn','-r0',strcat(SaveFilePath,'AllCellDiameter3DHistogram.tif'));
% % print(FigHandle01a,'-dtiffn','-r300',strcat(SaveFilePath,'AllCellDiameter3DHistogram.tif'));

% FigHandle01b=figure('Position',[40 300 600 600]);
% bar(CellDiameter3D_um_HistCtBinLoc,cumsum(CellDiameter3D_um_HistCt)./sum(CellDiameter3D_um_HistCt(:)),'FaceColor',[1 0.75 0],'BarWidth',1);
% ylabel(strcat('Cumulative Histogram: 3D Cell Diameter (um)'));
% xlabel('3D Cell Diameter (um)');
% xticks(CellDiameter3D_um_HistCtBinLoc);
% xtickangle(45);
% ax=gca;
% ax.FontSize=18;
% title({strcat('3D Cell Diameter (um)'),...
%     strcat('Total Cell Count (3D)=',num2str(size(CellDiameter3D_um_Input,1))),...
%     strcat('Radius/Sigma for Gaussian=',num2str(DiameterSigmaRatio_Query)),...
%     strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
% print(FigHandle01b,'-dtiffn','-r0',strcat(SaveFilePath,'AllCellDiameter3DCumHistogram.tif'));
% % print(FigHandle01b,'-dtiffn','-r300',strcat(SaveFilePath,'AllCellDiameter3DCumHistogram.tif'));
% close(FigHandle01b);

FigHandle01c=figure('Position',[50 300 600 600]);
% histogram(CellDiameter3D_um_Input,CellDiameter3D_um_HistEdge,'FaceColor',[1 0.75 0]);
bar(CellDiameter3D_um_HistCtBinLoc(MinIdx:end),CellDiameter3D_um_HistCt(MinIdx:end),'FaceColor',[0.75 1.00 0],'BarWidth',1);
ylabel('Histogram: 3D Cell Diameter (um)');
xlabel('3D Cell Diameter (um)');
xticks(CellDiameter3D_um_HistCtBinLoc(MinIdx:end));
xtickangle(45);
ax=gca;
ax.FontSize=18;
title({strcat('3D Cell Diameter (um)'),...
    strcat('Post Cutoff Cell Count (3D)=',num2str(sum(CellDiameter3D_um_HistCt(MinIdx:end)))),...
    strcat('Radius/Sigma for Gaussian=',num2str(DiameterSigmaRatio_Query)),...
    strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
% print(FigHandle01c,'-dtiffn','-r0',strcat(SaveFilePath,'PostCutCellDiameter3DHistogram.tif'));
print(FigHandle01c,'-dtiffn','-r150',strcat(SaveFilePath,'PostCutCellDiameter3DHistogram.tif'));

% FigHandle01d=figure('Position',[60 300 600 600]);
% bar(CellDiameter3D_um_HistCtBinLoc(MinIdx:end),cumsum(CellDiameter3D_um_HistCt(MinIdx:end))./sum(CellDiameter3D_um_HistCt(MinIdx:end)),'FaceColor',[0.75 1.00 0],'BarWidth',1);
% ylabel(strcat('Cumulative Histogram: 3D Cell Diameter>',num2str(MinCellDiameter_um,'%0.1f'),'(um)'));
% xlabel('3D Cell Diameter (um)');
% xticks(CellDiameter3D_um_HistCtBinLoc(MinIdx:end));
% xtickangle(45);
% ax=gca;
% ax.FontSize=18;
% title({strcat('3D Cell Diameter (um)'),...
%     strcat('Post Cutoff Cell Count (3D)=',num2str(sum(CellDiameter3D_um_HistCt(MinIdx:end)))),...
%     strcat('Radius/Sigma for Gaussian=',num2str(DiameterSigmaRatio_Query)),...
%     strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
% print(FigHandle01d,'-dtiffn','-r0',strcat(SaveFilePath,'PostCutCellDiameter3DCumHistogram.tif'));
% % print(FigHandle01d,'-dtiffn','-r300',strcat(SaveFilePath,'PostCutCellDiameter3DCumHistogram.tif'));
% close(FigHandle01d);

clearvars CellDiameter3D_um_Step CellDiameter3D_um_HistEdge_min CellDiameter3D_um_HistEdge_max;
clearvars FigHandle01 FigHandle01b;
clearvars Diff MinIdx;
clearvars *_Input;


%% Create Truncate CellCentroid3D_Diameter_um_LinIdx_Mtx based on MinCellDiameter_um

CellCentroid3D_Diameter_um_LinIdx_Mtx=CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx;
CellCentroid3D_Diameter_um_LinIdx_Mtx(CellCentroid3D_Diameter_um_LinIdx_Mtx(:,1)<MinCellDiameter_um,:)=[];
CellCentroid3D_Ct=size(CellCentroid3D_Diameter_um_LinIdx_Mtx,1);

clearvars CellCentroid3D_Diameter_um_LinIdx_AllCells_Mtx;


%% Prepare 3D Cell Position (Centroid) Map

fprintf('Preparing 3D Cell Position (Centroid) Map...\n');

CellCentroid3DStack=zeros(Img_Height,Img_Width,NumImgSlicesZSAR);
CellCentroid3DStack(CellCentroid3D_Diameter_um_LinIdx_Mtx(:,2))=1;

% imgaussfilt3 for display use only, sigma value not meaningful
CellCentroid3DStack=imgaussfilt3(CellCentroid3DStack,1.5);
CellCentroid3DStack=CellCentroid3DStack./max(CellCentroid3DStack(:))*1; % *2 to oversaturate
CellCentroid3DStack(CellCentroid3DStack>1)=1; % Intensity [0 1]


%% Display, Save 3D Cell Position (Centroid) Map

fprintf('Saving 3D Cell Position (Centroid) Map...\n');

SaveCellCentroid3DFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/CellCentroid3D_ZSAR/');
% mkdir(SaveCellCentroid3DFilePath);

% figure;
% for k=1:size(CellCentroid3DStack,3)
%     imshow(CellCentroid3DStack(:,:,k));
%     imwrite(double(CellCentroid3DStack(:,:,k)),strcat(SaveCellCentroid3DFilePath,'CellCentroid3D -',num2str(k,'%04.0f'),'.tif'));
% end

clearvars SaveCellCentroid3DFilePath;
clearvars k;


%% Overlay Cell Centroid (3D) Map on Original Image Stack

if OrgImgColocal_Query=='y'
    
    fprintf('Overlay 3D Cell Centroid Map on Original Image Stack...\n');
    
    CellCentroid3DStack_Input=CellCentroid3DStack;
    
    % Load orginal image stack for colocalization, if needed
    if OrgImgColocal_Query=='y'
        [OrgImgStack,OrgImg_Height,OrgImg_Width,NumOrgImgSlices]=ImgStackLoad_f(InputOrgImg_TifvsMAT_Query,BatchOrgImgInputFolder,OrgImgStackData_FilenameString,RGBChannel);
        if (NumOrgImgSlices~=NumImgSlices) || (OrgImg_Height~=Img_Height) || (OrgImg_Width~=Img_Width)
            error('Error. OrgImgStack size must equal ImgStack (GdMap) size.')
        end
    else
        OrgImgStack=[];
    end
    
    % Z-Shift, 1:1:1 Aspect Ratio Fix
    OrgImgStack=single(OrgImgStack(:,:,StartImgSliceZS:EndImgSliceZS)); % Z-Shift
    OrgImgStack=imresize3(single(OrgImgStack),[Img_Height Img_Width NumImgSlicesZSAR]); % Z-Shift, 1:1:1 Aspect Ratio
    OrgImgStack=OrgImgStack./max(OrgImgStack(:)); % Intensity [0 1]
    
    % Prepare overlay images and save
    SaveCellCentroid3DOrgImgColocalFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/CellCentroid3D_OrgImgStack_Colocal_ZSAR/');
    mkdir(SaveCellCentroid3DOrgImgColocalFilePath);
    
    CellCentroid3D_OrgImg_Colocal_Cell=cell(size(CellCentroid3DStack_Input,3),1);
    CellCentroid3D_OrgImg_Colocal_Cell(:)={zeros(Img_Height,Img_Width,3)};
    
    figure;
    for k=1:size(CellCentroid3DStack_Input,3)
        CellCentroid3D_OrgImg_Colocal_Cell{k,1}(:,:,1)=OrgImgStack(:,:,k)*1.50;   % Allow some over-saturation
        CellCentroid3D_OrgImg_Colocal_Cell{k,1}(:,:,2)=CellCentroid3DStack_Input(:,:,k);
        imshow(CellCentroid3D_OrgImg_Colocal_Cell{k,1});
        imwrite(double(CellCentroid3D_OrgImg_Colocal_Cell{k,1}),strcat(SaveCellCentroid3DOrgImgColocalFilePath,'CellCentroid3DOrgImgColocal -',num2str(k,'%04.0f'),'.tif'));
    end
    
end

clearvars OrgImgStack;
clearvars OrgImg_Height OrgImg_Width NumOrgImgSlices;
clearvars CellCentroid3D_OrgImg_Colocal_Cell;
clearvars CellCentroid3DStack;
clearvars SaveCellCentroid3DOrgImgColocalFilePath;
clearvars k;


%% Preparing Bone and Vessel Distance Stack: Load GdMap as Cell

% Load image stack (GdMap)
% All 3 colors
[ImgStack_R,~,~,~]=ImgStackLoad_20170328(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,1);
[ImgStack_G,~,~,~]=ImgStackLoad_20170328(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,2);
[ImgStack_B,~,~,~]=ImgStackLoad_20170328(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,3);

ImgStack_RGBCell=cell(3,1);
ImgStack_RGBCell{1,1}=single(ImgStack_R);
ImgStack_RGBCell{2,1}=single(ImgStack_G);
ImgStack_RGBCell{3,1}=single(ImgStack_B);

clearvars ImgStack_R ImgStack_G ImgStack_B;


%% Prepare bone, vessel and marrow cavity maps (with Z-Shift and 1:1:1 Aspect Ratio Fix)

fprintf('Preparing bone, vessel and marrow cavity maps, determine volumes...\n');

ImgStack_Cell_Input=ImgStack_RGBCell;

% == Bone Mask Volume
BoneLinIdx=intersect(intersect(find(abs(ImgStack_Cell_Input{1,1}-BoneColor(1,1))<0.5/255),...
    find(abs(ImgStack_Cell_Input{2,1}-BoneColor(1,2))<0.5/255)),...
    find(abs(ImgStack_Cell_Input{3,1}-BoneColor(1,3))<0.5/255));

BoneMaskStack=zeros(Img_Height,Img_Width,NumImgSlices,'logical');
BoneMaskStack(BoneLinIdx)=1;

BoneMaskStack=logical(BoneMaskStack(:,:,StartImgSliceZS:EndImgSliceZS)); % Z-Shift
BoneMaskStack=logical(imresize3(single(BoneMaskStack),[Img_Height Img_Width NumImgSlicesZSAR])); % Z-Shift, 1:1:1 Aspect Ratio

BoneVolume3D_px=numel(find(BoneMaskStack));
BoneVolume3D_mm=BoneVolume3D_px*((XY_PxLength.^2)*(Z_PxLength_ZSAR))*1E-9;

% == Vessel Mask Volume
VesselLinIdx=intersect(intersect(find(abs(ImgStack_Cell_Input{1,1}-VesselColor(1,1))<0.5/255),...
    find(abs(ImgStack_Cell_Input{2,1}-VesselColor(1,2))<0.5/255)),...
    find(abs(ImgStack_Cell_Input{3,1}-VesselColor(1,3))<0.5/255));

VesselMaskStack=zeros(Img_Height,Img_Width,NumImgSlices,'logical');
VesselMaskStack(VesselLinIdx)=1;

VesselMaskStack=logical(VesselMaskStack(:,:,StartImgSliceZS:EndImgSliceZS)); % Z-Shift
VesselMaskStack=logical(imresize3(single(VesselMaskStack),[Img_Height Img_Width NumImgSlicesZSAR])); % Z-Shift, 1:1:1 Aspect Ratio

VesselVolume3D_px=numel(find(VesselMaskStack));
VesselVolume3D_mm=VesselVolume3D_px*((XY_PxLength.^2)*(Z_PxLength_ZSAR))*1E-9;

% == BC Mask Volume
BCMaskLinIdx=intersect(intersect(find(abs(ImgStack_Cell_Input{1,1}-BCMaskColor(1,1))<0.5/255),...
    find(abs(ImgStack_Cell_Input{2,1}-BCMaskColor(1,2))<0.5/255)),...
    find(abs(ImgStack_Cell_Input{3,1}-BCMaskColor(1,3))<0.5/255));

BCMaskStack=ones(Img_Height,Img_Width,NumImgSlices,'logical');
BCMaskStack(BCMaskLinIdx)=0;  % Usual convention: 1 for Signal Areas

BCMaskStack=logical(BCMaskStack(:,:,StartImgSliceZS:EndImgSliceZS)); % Z-Shift
BCMaskStack=logical(imresize3(single(BCMaskStack),[Img_Height Img_Width NumImgSlicesZSAR])); % Z-Shift, 1:1:1 Aspect Ratio

BCMaskVolume3D_px=numel(find(~BCMaskStack));  % Volume of no-signal areas
BCMaskVolume3D_mm=BCMaskVolume3D_px*((XY_PxLength.^2)*(Z_PxLength_ZSAR))*1E-9;

% == Img Mask Volume
ImgMaskLinIdx=intersect(intersect(find(abs(ImgStack_Cell_Input{1,1}-ImgMaskColor(1,1))<0.5/255),...
    find(abs(ImgStack_Cell_Input{2,1}-ImgMaskColor(1,2))<0.5/255)),...
    find(abs(ImgStack_Cell_Input{3,1}-ImgMaskColor(1,3))<0.5/255));

ImgMaskStack=ones(Img_Height,Img_Width,NumImgSlices,'logical');
ImgMaskStack(ImgMaskLinIdx)=0;  % Usual convention: 1 for Signal Areas

ImgMaskStack=logical(ImgMaskStack(:,:,StartImgSliceZS:EndImgSliceZS)); % Z-Shift
ImgMaskStack=logical(imresize3(single(ImgMaskStack),[Img_Height Img_Width NumImgSlicesZSAR])); % Z-Shift, 1:1:1 Aspect Ratio

ImgMaskVolume3D_px=numel(find(~ImgMaskStack));  % Volume of no-signal areas
ImgMaskVolume3D_mm=ImgMaskVolume3D_px*((XY_PxLength.^2)*(Z_PxLength_ZSAR))*1E-9;

% == Marrow Cavity Volume
MarrowCavityVolume3D_px=numel(ImgMaskStack)-BoneVolume3D_px-VesselVolume3D_px-ImgMaskVolume3D_px-BCMaskVolume3D_px;
MarrowCavityVolume3D_mm=MarrowCavityVolume3D_px*((XY_PxLength.^2)*(Z_PxLength_ZSAR))*1E-9;


clearvars BoneLinIdx VesselLinIdx BCMaskLinIdx ImgMaskLinIdx;
clearvars ImgStack_RGBCell;
clearvars BoneVolume3D_px VesselVolume3D_px ImgMaskVolume3D_px MarrowCavityVolume3D_px;
clearvars *_Input;


%% 3D Euclidean Distance for Bone and Vessel Mask Maps

% fprintf('Calculate 3D Euclidean Distance for bone and vessel mask stacks...\n');
fprintf('Calculate 3D Euclidean Distance for  vessel mask stack...\n');

BoneMaskStack_Input=BoneMaskStack;
VesselMaskStack_Input=VesselMaskStack;

% Calculate 3D euclidean distance to bone and vessel
% BoneDistLinIdxStack records linear index of closest bone for each
% coordinate in stack
[BoneDistMapStack,BoneDistLinIdxStack]=bwdist(BoneMaskStack_Input,'euclidean');
[VesselDistMapStack,VesselDistLinIdxStack]=bwdist(VesselMaskStack_Input,'euclidean');

clearvars BoneMaskStack VesselMaskStack;
clearvars *_Input;


%% 3D cell distance to bone and vasculature
% Inputs must have 1:1:1 voxel XYZ aspect ratio for accurate measurement

% fprintf('Plot: Distance from 3D Cells to Bone and Vasculature...\n');
fprintf('Plot: Distance from 3D Cells to Vasculature...\n');

CellCentroid3D_LinIdx_Input=CellCentroid3D_Diameter_um_LinIdx_Mtx(:,2);
BoneDistMapStack_Input=BoneDistMapStack;
VesselDistMapStack_Input=VesselDistMapStack;
BoneDistLinIdxStack_Input=BoneDistLinIdxStack;
VesselDistLinIdxStack_Input=VesselDistLinIdxStack;

% = Include option for all non-bone/vessel pixels, or random pixels
% Do not allow less than 1 cell radius from bone or vessel surface
BoneDistMapStack_MinCellSizeExclude=BoneDistMapStack_Input;
VesselDistMapStack_MinCellSizeExclude=VesselDistMapStack_Input;
% BoneDistMapStack_MinCellSizeExclude(BoneDistMapStack_MinCellSizeExclude<MinCellDiameter_um/XY_PxLength*0.5)=0;
% VesselDistMapStack_MinCellSizeExclude(VesselDistMapStack_MinCellSizeExclude<MinCellDiameter_um/XY_PxLength*0.5)=0;
AllCentroid3D_LinIdx=...
    intersect(...
    intersect(...
    intersect(find(BoneDistMapStack_MinCellSizeExclude),find(VesselDistMapStack_MinCellSizeExclude)),...
    find(BCMaskStack)),...
    find(ImgMaskStack));
RdmCentroid3D_LinIdx_Pick=round(rand(size(CellCentroid3D_LinIdx_Input,1),1)*size(AllCentroid3D_LinIdx,1));
RdmCentroid3D_LinIdx=AllCentroid3D_LinIdx(RdmCentroid3D_LinIdx_Pick,1);

% = Create cell to hold distance data
% Cell row for CellPos type:
% Row 1: CellCentroid3D
% Row 2: All cavity pixels
% Row 3: Random cavity pixels
% Cell col for distance and histogram data:
% Col 1: Dist3D_um
% Col 2: HistCt
% Col 2: HistEdge

BoneDist3DData_um_Cell=cell(3,3);
VesselDist3DData_um_Cell=cell(3,3);

for CellPosSet=1:1
    
    if CellPosSet==1    % Cell Centroid 3D
        CellCentroid3D_LinIdx_Input=CellCentroid3D_Diameter_um_LinIdx_Mtx(:,2);
        HistIDString=',  CellCentroid3D';
        SaveDistHistFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/CellCentroid3D_DistHistgm/');
    elseif CellPosSet==2  % All non bone/vessel/mask cavity pixels
        CellCentroid3D_LinIdx_Input=AllCentroid3D_LinIdx;
        HistIDString=',  AllCavPixels';
        SaveDistHistFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/AllCavPx_DistHistgm/');
    elseif CellPosSet==3  % Random non bone/vessel/mask cavity pixels
        CellCentroid3D_LinIdx_Input=RdmCentroid3D_LinIdx;
        HistIDString=',  RdmCavPixels';
        SaveDistHistFilePath=strcat(BatchImgOutputFolder,'_',timestamp,'/RdmCavPx_DistHistgm/');
    end
    
    mkdir(SaveDistHistFilePath);
    
    % Get distance to bone or vessel for each 3D cell, convert to um
    CellBoneDist3D_um=BoneDistMapStack_Input(CellCentroid3D_LinIdx_Input)*...
        ((XY_PxLength*2+Z_PxLength_ZSAR)/3);
    CellVesselDist3D_um=VesselDistMapStack_Input(CellCentroid3D_LinIdx_Input)*...
        ((XY_PxLength*2+Z_PxLength_ZSAR)/3);
    
    % Get 99.9% percentile, for limiting x axis during plotting
    CellBoneDist3D_um_xMax=prctile(CellBoneDist3D_um,99.9);
    CellVesselDist3D_um_xMax=prctile(CellVesselDist3D_um,99.9);
    
    % Get LinIdx of bone or vessel where distance is recorded for each cell
    CellBoneDist3D_LinIdx=BoneDistLinIdxStack_Input(CellCentroid3D_LinIdx_Input);
    CellVesselDist3D_LinIdx=VesselDistLinIdxStack_Input(CellCentroid3D_LinIdx_Input);
    
    Dist3D_um_step=2.5;
    
    BoneDist3D_um_HistEdge=0:Dist3D_um_step:ceil(max(CellBoneDist3D_um./Dist3D_um_step))*Dist3D_um_step;
    [BoneDist3D_um_HistCt,~]=histcounts(CellBoneDist3D_um,BoneDist3D_um_HistEdge);  % For cumulative histogram, saving
    BoneDist3D_um_HistCtBinLoc=0.5*(BoneDist3D_um_HistEdge(1:end-1)+BoneDist3D_um_HistEdge(2:end)); % For cumulative histogram
    
    VesselDist3D_um_HistEdge=0:Dist3D_um_step:ceil(max(CellVesselDist3D_um./Dist3D_um_step))*Dist3D_um_step;
    [VesselDist3D_um_HistCt,~]=histcounts(CellVesselDist3D_um,VesselDist3D_um_HistEdge);
    VesselDist3D_um_HistCtBinLoc=0.5*(VesselDist3D_um_HistEdge(1:end-1)+VesselDist3D_um_HistEdge(2:end)); % For cumulative histogram
    
    BoneDist3DData_um_Cell{CellPosSet,1}=CellBoneDist3D_um;
    BoneDist3DData_um_Cell{CellPosSet,2}=BoneDist3D_um_HistCt;
    BoneDist3DData_um_Cell{CellPosSet,3}=BoneDist3D_um_HistEdge;
    
    VesselDist3DData_um_Cell{CellPosSet,1}=CellVesselDist3D_um;
    VesselDist3DData_um_Cell{CellPosSet,2}=VesselDist3D_um_HistCt;
    VesselDist3DData_um_Cell{CellPosSet,3}=VesselDist3D_um_HistEdge;
    
% % %     FigHandle02=figure('Position',[110 300 600 600]);
% % %     bar(BoneDist3D_um_HistCtBinLoc,100*BoneDist3D_um_HistCt./sum(BoneDist3D_um_HistCt(:)),'FaceColor',[0 0 1],'BarWidth',1);
% % %     xlim([0 CellBoneDist3D_um_xMax]);
% % %     ylabel(strcat('Histogram: Cell Ct %',HistIDString));
% % %     xlabel('3D Distance to Nearest Bone (um)');
% % %     ytickformat('percentage');
% % %     ax=gca;
% % %     ax.FontSize=18;
% % %     title({strcat('3D Cell Distance to Bone (um)'),...
% % %         strcat('Post cutoff cell count (3D)=',num2str(size(CellBoneDist3D_um,1))),...
% % %         strcat('Bone Volume (mm3)=',num2str(BoneVolume3D_mm,'%0.4E')),...
% % %         strcat('Marrow Cavity Volume (mm3)=',num2str(MarrowCavityVolume3D_mm,'%0.4E')),...
% % %         strcat('Bone Volume Fraction=',num2str(BoneVolume3D_mm/MarrowCavityVolume3D_mm,'%0.4E')),...
% % %         strcat('Cell Density, Marrow Cavity (#/mm3)=',num2str(size(CellBoneDist3D_um,1)/MarrowCavityVolume3D_mm,'%0.4E')),...
% % %         strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
% % %     print(FigHandle02,'-dtiffn','-r300',strcat(SaveDistHistFilePath,'CellBone3DDistHistogram.tif'));
% % %     if CellPosSet>1
% % %         close(FigHandle02);
% % %     end
    
    FigHandle03=figure('Position',[160 300 600 600]);
    bar(VesselDist3D_um_HistCtBinLoc,100*VesselDist3D_um_HistCt./sum(VesselDist3D_um_HistCt(:)),'FaceColor',[1 0 0],'BarWidth',1);
    xlim([0 CellVesselDist3D_um_xMax+(VesselDist3D_um_HistCtBinLoc(2)-VesselDist3D_um_HistCtBinLoc(1))]);
    ylabel(strcat('Histogram: Cell Ct %',HistIDString));
    xlabel('3D Distance to Nearest Vessel (um)');
    ytickformat('percentage');
    ax=gca;
    ax.FontSize=18;
    title({strcat('3D Cell Distance to Vasculature (um)'),...
        strcat('Post cutoff cell count (3D)=',num2str(size(CellVesselDist3D_um,1))),...
        strcat('Vessel Volume (mm3)=',num2str(VesselVolume3D_mm,'%0.4E')),...
        strcat('Marrow Cavity Volume (mm3)=',num2str(MarrowCavityVolume3D_mm,'%0.4E')),...
        strcat('Vessel Volume Fraction=',num2str(VesselVolume3D_mm/MarrowCavityVolume3D_mm,'%0.4E')),...
        strcat('Cell Density, Marrow Cavity (#/mm3)=',num2str(size(CellVesselDist3D_um,1)/MarrowCavityVolume3D_mm,'%0.4E')),...
        strcat('MinCellDiameter-um-Query=',num2str(MinCellDiameter_um,'%0.1f'))},'FontSize',11);
    % print(FigHandle03,'-dtiffn','-r0',strcat(SaveDistHistFilePath,'CellVessel3DDistHistogram.tif'));
    print(FigHandle03,'-dtiffn','-r150',strcat(SaveDistHistFilePath,'CellVessel3DDistHistogram.tif'));
    if CellPosSet>1
        close(FigHandle03);
    end
        
    % Save Distance Maps
% % %     SaveBoneDistMap3DFilePath=...
% % %         strcat(SaveDistHistFilePath,'/BoneDistMap3D 0-',num2str(max(BoneDistMapStack_Input(:))*((XY_PxLength*2+Z_PxLength_ZSAR)/3)),'um/');
% % %     mkdir(SaveBoneDistMap3DFilePath);
    
% % %     SaveVesselDistMap3DFilePath=...
% % %         strcat(SaveDistHistFilePath,'/VesselDistMap3D 0-',num2str(max(VesselDistMapStack_Input(:))*((XY_PxLength*2+Z_PxLength_ZSAR)/3)),'um/');
% % %     mkdir(SaveVesselDistMap3DFilePath);
    
% % %     SaveCellPosBoneDistMap3DFilePath=...
% % %         strcat(SaveDistHistFilePath,'/CellPosBoneDistMap3D 0-',num2str(max(BoneDistMapStack_Input(:))*((XY_PxLength*2+Z_PxLength_ZSAR)/3)),'um/');
% % %     mkdir(SaveCellPosBoneDistMap3DFilePath);
    
% % %     SaveCellPosVesselDistMap3DFilePath=...
% % %         strcat(SaveDistHistFilePath,'/CelPosVesselDistMap3D 0-',num2str(max(VesselDistMapStack_Input(:))*((XY_PxLength*2+Z_PxLength_ZSAR)/3)),'um/');
% % %     mkdir(SaveCellPosVesselDistMap3DFilePath);
    
    % figure;
% % %     for k=1:size(BoneDistMapStack_Input,3)
% % %         % imshow(double(BoneDistMapStack_Input(:,:,k)/max(BoneDistMapStack_Input(:))));
% % %         imwrite(double(BoneDistMapStack_Input(:,:,k)/max(BoneDistMapStack_Input(:))),strcat(SaveBoneDistMap3DFilePath,'BoneDistMap3D -',num2str(k,'%04.0f'),'.tif'));
% % %     end
    
    % figure;
% % %     for k=1:size(VesselDistMapStack_Input,3)
% % %         % imshow(double(VesselDistMapStack_Input(:,:,k)/max(VesselDistMapStack_Input(:))));
% % %         imwrite(double(VesselDistMapStack_Input(:,:,k)/max(VesselDistMapStack_Input(:))),strcat(SaveVesselDistMap3DFilePath,'VesselDistMap3D -',num2str(k,'%04.0f'),'.tif'));
% % %     end
    
    CellCentroid3D_LinIdx_Input_Stack=zeros(size(BoneDistMapStack_Input),'single');
    CellCentroid3D_LinIdx_Input_Stack(CellCentroid3D_LinIdx_Input)=1;
    CellCentroid3D_LinIdx_Input_Stack=imgaussfilt3(single(CellCentroid3D_LinIdx_Input_Stack),1.5);
    CellCentroid3D_LinIdx_Input_Stack=CellCentroid3D_LinIdx_Input_Stack./max(CellCentroid3D_LinIdx_Input_Stack(:));
    
% % %     % figure;
% % %     for k=1:size(BoneDistMapStack_Input,3)
% % %         Img=zeros(size(BoneDistMapStack_Input,1),size(BoneDistMapStack_Input,2),3);
% % %         Img(:,:,1)=double(BoneDistMapStack_Input(:,:,k)/max(BoneDistMapStack_Input(:)));
% % %         Img(:,:,2)=double(CellCentroid3D_LinIdx_Input_Stack(:,:,k));
% % %         % imshow(double(Img));
% % %         imwrite(double(Img),strcat(SaveCellPosBoneDistMap3DFilePath,'CellPosBoneDistMap3D -',num2str(k,'%04.0f'),'.tif'));
% % %     end
    
    % figure;
% % %     for k=1:size(VesselDistMapStack_Input,3)
% % %         Img=zeros(size(VesselDistMapStack_Input,1),size(VesselDistMapStack_Input,2),3);
% % %         Img(:,:,1)=double(VesselDistMapStack_Input(:,:,k)/max(VesselDistMapStack_Input(:)));
% % %         Img(:,:,2)=double(CellCentroid3D_LinIdx_Input_Stack(:,:,k));
% % %         % imshow(double(Img));
% % %         imwrite(double(Img),strcat(SaveCellPosVesselDistMap3DFilePath,'CellPosVesselDistMap3D -',num2str(k,'%04.0f'),'.tif'));
% % %     end
    
    % Save distance distributions
% % %     save(strcat(SaveDistHistFilePath,'BoneDist3Ddata.mat'),'BoneDist3D_um_HistCt','BoneDist3D_um_HistCtBinLoc');
    save(strcat(SaveDistHistFilePath,'VesselDist3Ddata.mat'),'VesselDist3D_um_HistCt','VesselDist3D_um_HistCtBinLoc');
    
end

clearvars AllCentroid3D_LinIdx RdmCentroid3D_LinIdx;
clearvars BoneDistMapStack VesselDistMapStack;
clearvars BoneDistMapStack_MinCellSizeExclude;
clearvars BoneDistLinIdxStack VesselDistLinIdxStack;
clearvars BoneDist3D_um_HistCtBinLoc VesselDist3D_um_HistCtBinLoc;
clearvars RdmCentroid3D_LinIdx_Pick;
clearvars CellBoneDist3D_um CellVesselDist3D_um;
clearvars BoneDist3D_um_HistCt VesselDist3D_um_HistCt;
clearvars BoneDist3D_um_HistEdge VesselDist3D_um_HistEdge;
clearvars BoneDistMapStack_MinCellSizeExclude VesselDistMapStack_MinCellSizeExclude;
clearvars ImgMaskStack BCMaskStack;
clearvars CellCentroid3D_LinIdx_Input_Stack;
clearvars Dist3D_um_step;
clearvars CellRadiusMin CellRadiusMax;
clearvars FigHandle*;
clearvars Img k;
clearvars *_Input;


%% Save script in directory
% html format; do not evaulate code or save figures

% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);


%% Save workspace

clearvars Batch*Folder;
clearvars *FileID;
clearvars *_FilenameString;
clearvars Save*FilePath -except SaveFilePath;
close all;
% save(strcat(SaveFilePath,'Workspace_3DCellCtEuclidDist.mat'));

