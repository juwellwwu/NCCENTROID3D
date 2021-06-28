function [ImgStack_XYAv]=ImgStackXYAv(ImgStack_Input,ImgMask_Query,ImgMaskStack,FlatStrel_Query,strel_Radius,offsetstrel_Radius)

fprintf('Running XY average filter...\n');

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

%  Apply mask on ImgStack for averaging
if ImgMask_Query=='y'
     ImgStack_Input(~ImgMaskStack)=NaN;
end

% Average filter, n-pixel radius circle
if FlatStrel_Query=='y'
    h_av=fspecial('disk',strel_Radius);
else
    h_av=fspecial('disk',offsetstrel_Radius);
end

% Average filter on each Z slice.
% imfilter() is okay with NaN; if neighbourhood touches NaN pixels, 
% Center pixel = NaN 
ImgStack_XYAv=single(zeros(size(ImgStack_Input)));
for k=1:NumImgSlices
    ImgStack_XYAv(:,:,k)=single(imfilter(ImgStack_Input(:,:,k),h_av));
end

% Set ImgMask areas to 0
if ImgMask_Query=='y'
    ImgStack_XYAv(~ImgMaskStack)=0;
end

end

