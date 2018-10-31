clear all 
close all
clc

path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

startim = 26;
endim = 38;
total_img_position = 5;
irt_img_basename = 'IRT_0';
rgb_img_basename = 'RGB_0';
leg = 'l';

plotIndex = 0;
for p=startim:endim
    if p>9
      irt_img_basename = 'IRT_';
      rgb_img_basename = 'RGB_'; 
    end
    
    for k=1:total_img_position
        
        if k==1
            leg = 'l';
        elseif k==2
            leg = 'r';
        elseif k==3
            leg = 'f';
        elseif k==4
            leg = 'i';
        else
            leg = 'b';
        end
        
    irt_img = [irt_img_basename,num2str(p),'_',leg]
    rgb_img = [rgb_img_basename,num2str(p),'_',leg]
    
    
    fn = fullfile(strcat(path,'\',[irt_img,'.asc']));
    
    if exist(fn,'file')==0
        % scan another file
        continue
    end
     

irt_width = 1024; 
irt_height = 768;



% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));

% Convert RGB image to chosen color space
I = rgb2lab(rgb_img);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 10.355;
channel1Max = 93.949;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 3.473;
channel2Max = 48.305;

% Define thresholds for channel 3 based on histogram settings
channel3Min = -30.183;
channel3Max = 35.075;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = rgb_img;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

BW = bwareaopen(BW,100);
BW = imfill(BW,'holes');
% Find edges using the Canny operator with hysteresis thresholds of 0.1
% and 0.2 with smoothing parameter sigma set to 1.
edgeim_rgb = edge(BW,'canny', [0.1 0.9], 10);


plotIndex = plotIndex +1;
figure
% subplot(2,2,plotIndex)
imshow(edgeim_rgb);
    end
end