clear all 
close all
clc

path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img_basename = 'IRT_0';
rgb_img_basename = 'RGB_0';
leg = '_l';
count = 0;


for p=[33,36,38]    % scan the image patient one by one
    if p>9
      irt_img_basename = 'IRT_';
      rgb_img_basename = 'RGB_'; 
    end
    
    for k=1:5
        
        if k==1
            leg = '_l';
        elseif k==2
            leg = '_r';
        elseif k==3
            leg = '_f';
        elseif k==4
            leg = '_i';
        else
            leg = '_b';
        end
        
    irt_img = [irt_img_basename,num2str(p),leg]
    rgb_img = [rgb_img_basename,num2str(p),leg]
    
    
    fn = fullfile(strcat(path,'\',[irt_img_basename,num2str(p),leg,'.asc']));
    
    % check whether the image exist
    if exist(fn,'file')==0
        % if not scan another file
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

% % Compress/rescale the irt image
% [rows, columns, numberOfColorChannels] = size(irt_img);
% newWidth = round(0.87 * columns)
% stretchedImage = imresize(irt_img, [rows newWidth]);
% irt_img = stretchedImage;

% Uncomment this if want to select the movingPoints and fixedPoints
leg_fixed = rgb_img;
leg_moving = irt_img; 
f = cpselect(leg_moving, leg_fixed);
uiwait(f);

% fixedPoints = [6.031250000000001e+02 2.461249999999999e+02;8.903750000000002e+02 1.103749999999999e+02;6.128750000000002e+02 6.181250000000000e+02;7.688750000000001e+02 6.361250000000000e+02];
% movingPoints = [1.531250000000000e+02 2.078749999999999e+02;4.328750000000001e+02 22.624999999999886;2.348750000000001e+02 7.343750000000000e+02;4.021250000000001e+02 7.501250000000000e+02];

count = count +1;

movingPointsArray{1,count} = movingPoints;
fixedPointsArray{1,count} = fixedPoints;

% to ensure no duplicate variable 
clear movingPoints;
clear fixedPoints;

% Create a geometric transformation that can be used to align the two images
tform{1,count} = fitgeotrans(movingPointsArray{1,count},fixedPointsArray{1,count},'affine');


    end
end
 

%% find mean of each transformation matix
tform_mean=0;
for c=1:count
    tform_mean = tform_mean+tform{1,c}.T;
end
tform_mean = tform_mean./count;
%tform_mean = affine2d(tform_mean);