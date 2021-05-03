function [out, ratio] = ratioChan_LCv2_mask_comp_DCMGlut(img)


Show_Channels  = 1;    % If 0 you will see less images through the process

Edge_Draw      = 0;    % Select this if you want to draw in the cells edge for distance to edge calculation

convf          = 0.09; % micron to pixel conversion factor

Size_Low       = 20;   % Low end size constraint

Size_High      = 1e5;  % High end size constraint

Channel2_scale = 1;    % Channel 2 scale factor

Channel1_scale = 1;  % Channel 1 scale factor

Thresh_combine = 0;    % 0: Threshold both images separately, then combine the mask. 1: combine the images then threshold

Open_Close     = 0;    % Watershed on/off

pwr            = 3;    % Watershed Strength

ColormapScale  = 2;    % Sets the colormap max value; original value 8. Too high?

ColormapBins   = 20;   % This sets the number of bins in the colormap

MaxPixel       = 500;  % This is the maximum intensity of the pixels in image default 65000
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
close all
%edited to load custom colormap
% colormap(hot(ColormapBins));
% currentMap = colormap;
% newMap = [0 0 0; currentMap];
newMap=inferno(ColormapBins); %only new addition

%InitializeImg = imread(img);
%Reading in the images:
img1 = mat2gray(imread(img,2),[0 MaxPixel]); %Cy5 (H5) channel - scale to pixel intensity
img2 = mat2gray(imread(img,1),[0 MaxPixel]); %RFP (Fn) channel - scale to pixel intensity

% account for third matrix dimension
X = imread(img);
[~,~,h] = size(X);
if h == 3
    img1 = sum(img1,3)./3;
    img2 = sum(img2,3)./3;
end

compimg = img1+img2;

figure;
if Show_Channels == 1;
    
subplot(1,3,1)
colormap(newMap)
imagesc(img1);
title('H5 Channel')

subplot(1,3,2)
colormap(newMap)
imagesc(img2);
title('Fn Channel')

end

%% A binary mask, using Otsu thresholding 

if Thresh_combine == 0
    
    threshold1 = graythresh(img1); %This is the Otsu threshold for img1
    mask1 = im2bw(img1,threshold1); %This mask is constructed of the first chanel
    threshold2 = graythresh(img2); %This is the Otsu threshold for img2
    mask2 = im2bw(img2,threshold2); %This mask is constructed of the second channel
    
    if Open_Close == 1;
        
        mask1 = OpenCloseImage(mask1,pwr);
        
        mask1 = double(mask1);
        
        mask2 = OpenCloseImage(mask2,pwr);
        
        mask2 = double(mask2);
        
    end
    
    figure;
    colormap(newMap);
    mask = mask1~= false & mask2 ~= false; % mask to only accept pixels that are nonzero in both image threshold masks
    imagesc(mask);
    
elseif Thresh_combine == 1
    
    threshold = graythresh(compimg); %This is the Otsu threshold for both channels combined
    mask = im2bw(compimg,threshold); %This mask is constructed of a composite of both channels
    
    if Open_Close == 1;
        
        mask = OpenCloseImage(mask,pwr);
        
        mask = double(mask);
        
    end
    
    figure;
    colormap(newMap);
    imagesc(mask);
    
else
    error('Set Thresh_combine to 0 or 1');
end
%% Analysis on a pixel by pixel basis
img22 = img2.*mask;
img11 = img1.*mask;

img22 = img22.*Channel2_scale;      % You can scale channels here
img11 = img11.*Channel1_scale;

figure
subplot(1,2,1)
img22(isnan(img22)) = 0;
colormap(newMap);
imagesc(img22);
title('Channel Two w/ Mask')

subplot(1,2,2)
img11(isnan(img11))=0;
colormap(newMap);
imagesc(img11);
title('Channel One w/ Mask')

%% Now lets choose an ROI for further analysis
ROIselection = roipoly;

img22 = img2.*mask;
img11 = img1.*mask;

img22(~ROIselection) = 0;
img11(~ROIselection) = 0;

img22 = img22.*Channel2_scale;      % You can scale channels here
img11 = img11.*Channel1_scale;

%subplot(1,2,1)
%img22(isnan(img22)) = 0;
%colormap(newMap);
%imagesc(img22);
%title('Channel Two w/ ROI and Mask')

subplot(1,1,1)
img11(isnan(img11))=0;
colormap(newMap);
imagesc(img11);
title('Channel One w/ ROI and Mask')
%% Now we can Ratio the two channels and visualize

rat = img11./img22;
rat(isnan(rat)) = 0;
figure;
colormap(newMap);
imagesc(rat); 
caxis([0 ColormapScale]);
% title('Ratio of Channel H5/Fn (Selections)')
colorbar('AxisLocation','out')
axis equal
axis off

ratb = rat;
ratb(ratb>0) = 1;
ratb(ratb~=1) = 0;

%% Finding the overall median and mean value of the ratioed chanels

poi = rat(rat ~= 0);
ratio = poi;
medrat = median(poi);
avgrat = mean(poi);
display(medrat,'Cy5:RFP median')
display(avgrat,'Cy5:RFP average')


%% New lets capture some parameters about the ratioed channels

S = bwconncomp(ratb);
data = regionprops(S,'Area');
  
% Area Threshold 
idx = find([data.Area] > Size_Low & [data.Area]<=Size_High); 
aimg = ismember(labelmatrix(S),idx);

S2 = bwconncomp(mat2gray(aimg));
data = regionprops(S2,rat,'Area','MajorAxis','MinorAxis','Orientation','Solidity','MeanIntensity', ...
                   'MaxIntensity');            

% calculate the median values of each pixel and add to structure array
for a = 1:length(S2.PixelIdxList) % run through each ROI
    indices = (S2.PixelIdxList{a}); 
    vals = rat([indices]);
    data(a).MedianIntensity = median(vals);   
end

%% Multiplication of Conversion Factior (Convf) 

for k = 1:numel(data);
      data(k).Area = data(k).Area.*convf.^2;
      data(k).MinorAxisLength = data(k).MinorAxisLength.*convf;
      data(k).MajorAxisLength = data(k).MajorAxisLength.*convf;
end
               
out = data;

