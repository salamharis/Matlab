clear all 
close all
clc

path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img_basename = 'IRT_0';
rgb_img_basename = 'RGB_0';
leg = 'l';

for p=1:5
    if p>9
      irt_img_basename = 'IRT_';
      rgb_img_basename = 'RGB_'; 
    end
    
    for k=1:5
        
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
     

% get 2D matrix (irt_img) containing temperature values for each pixel in °C
irt_img = dir(strcat(path,'\',irt_img,'*.asc'));
irt_img = textscan(fopen(irt_img(1,1).name),'%q');

irt_width = irt_img{1,1}{3,1}; 
irt_width = str2double(irt_width(12:length(irt_width))); 
irt_height = irt_img{1,1}{4,1};
irt_height = str2double(irt_height(13:length(irt_height)));

image = zeros(irt_height, irt_width);
for j=16:size(irt_img{1,1},1)
    w = 1 + mod((j-16),irt_width);
    h = floor((j-16)/irt_width + 1);
    image(h,w) = str2double(strrep(irt_img{1,1}{j,1},',','.'));
end
irt_img = image;
irt_img = (irt_img - min(irt_img(:)))/(max(irt_img(:)) - min(irt_img(:)));


% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));

% Compress/rescale the irt image
[rows, columns, numberOfColorChannels] = size(irt_img);
newWidth = round(0.87 * columns)
stretchedImage = imresize(irt_img, [rows newWidth]);
irt_img = stretchedImage;

% treshold the irt img to highlight the high gradient region
b_region = zeros(irt_height, irt_width);
threshold_value = 0.5;
for w=1:newWidth
    for h=1:irt_height
        if irt_img(h,w) <= threshold_value
            b_region(h,w) = 0;
        else
            b_region(h,w) = irt_img(h,w);
        end
    end
end
irt_img = b_region;

% insert the transformation function
% 1 till 10
% tform_mean = [0.855948865950106 0.0973273435355614 0;-0.0973273435355614 0.855948865950106 0;436.986969263621 -11.6329975478634 1];
%10 till 20
% tform_mean = [0.826759554004235 0.0527430881581958 0;-0.0527430881581958 0.826759554004235 0;401.684329628535 17.1661710177828 1];
tform_mean = [0.873508357250477    0.051630954289839   0; ...
             -0.101630954289839    0.853508357250477   0; ...
             4.555054787023003e+02 -15.164797761514200 1];
         
tform_mean = affine2d(tform_mean);

irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
img_title = ['Image ',num2str(p),' position ',leg];
figure
imshowpair(rgb_img,irt_registered,'blend')
title(img_title);
    end
end
