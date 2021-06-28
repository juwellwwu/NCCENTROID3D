function [ImgStack_XYMin]=ImgStackXYMin(ImgStack_Input,ImgMask_Query,ImgMaskStack,FlatStrel_Query,strel_Radius,offsetstrel_Radius,offsetstrel_MaxOffset)

fprintf('Running XY minimum filter...\n');

if isa(ImgStack_Input,'single')
    SatVal=1;
elseif isa(ImgStack_Input,'double')
    SatVal=1;
    ImgStack_Input=single(ImgStack_Input);
else
    error('Error. ImgStack_Input should be of class SINGLE or DOUBLE.')
end

Img_Height=size(ImgStack_Input,1);
Img_Width=size(ImgStack_Input,2);
NumImgSlices=size(ImgStack_Input,3);

% Prepare for min filtering
% strel creates flat structural element for neighbourhood
% offsetstrel creates non-flat structural element for neighbourhood
if FlatStrel_Query=='y'
    SE=strel('disk',strel_Radius,8);
else
    SE=offsetstrel('ball',offsetstrel_Radius,offsetstrel_MaxOffset,8);
end

% imerode does not take NaN; set NaN to 1 (largely does not affect min())
ImgStack_XYMin=single(ImgStack_Input);
ImgStack_XYMin(isnan(ImgStack_XYMin))= 1;

% Min filter on each Z slice.
% MATLAB's imerode is the same as min filter for grayscale images
parfor k=1:NumImgSlices
    ImgStack_XYMin(:,:,k)=single(imerode(ImgStack_XYMin(:,:,k),SE));
end

% Set ImgMask areas to 0
if ImgMask_Query=='y'
    ImgStack_XYMin(~ImgMaskStack)=0;
end

end
