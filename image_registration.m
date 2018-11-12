clear all 
close all
clc

path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

num_of_img = 21;
total_img_position = 5;
irt_img_basename = 'IRT_0';
rgb_img_basename = 'RGB_0';
leg = 'l';

for p=26:33
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

% treshold the irt img to highlight the high gradient region
% b_region = zeros(irt_height, irt_width);
% threshold_value = 0.5;
% for w=1:irt_width
%     for h=1:irt_height
%         if irt_img(h,w) <= threshold_value
%             b_region(h,w) = 0;
%         else
%             b_region(h,w) = irt_img(h,w);
%         end
%     end
% end
% irt_img = b_region;

% insert the transformation function
%     [ScaleX   ShearY   0;
%      ShearX   ScaleY   0;
%      TranslaX TranslaY 1]

% tform_mean = [0.773508357250477    0.051630954289839    0; ...
%              -0.081630954289839    0.803508357250477    0; ...
%              4.505054787023003e+02 25.164797761514200   1];

%4-21
% tform_mean = [0.805599525356956 0.020598075940177 0;-0.104893577450360 0.781929514272011 0;4.302242013663251e+02 48.320819733806360 1];
%25-27
% tform_mean = [0.637487839738744 -0.004434008430408 0;-0.040811608158753 0.551993295992885 0;4.416109116294455e+02 1.549737714257566e+02 1];
%28-32
% tform_mean = [0.878866245374861 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;3.605920167976121e+02 52.808982154224054 1];
%33-38
tform_mean =[0.822833383351964 8.876458614865464e-04 0;-0.023691893394736 0.829565435596974 0;2.270445581698194e+02 56.856631985064155 1];
tform_mean = affine2d(tform_mean);

irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
img_title = ['Image ',num2str(p),' position ',leg];
figure
imtool(irt_img);
% imshowpair(rgb_img,irt_registered,'blend')
title(img_title);
    end
end
