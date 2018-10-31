path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img = 'IRT_11_l';
rgb_img = 'RGB_11_l';

irt_img_bmp = imread(strcat(path,'\',irt_img,'.bmp'));


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

b_region = zeros(irt_height, irt_width);
threshold_value = 0.5;
for w=1:irt_width
    for h=1:irt_height
        if irt_img(h,w) <= threshold_value;
            b_region(h,w) = 0;
        else
            b_region(h,w) = irt_img(h,w);
        end
    end
end
figure
imshow(b_region)
% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));
irt_img_bmp = imresize(irt_img_bmp,min(irt_width,irt_height)/min(size(irt_img_bmp(:,:,1))));

% crop the right side parameter
% irt_img_bmp = imcrop(irt_img_bmp,[1, 1, 1023, 768]);
% imtool(irt_img_bmp)
% imtool(rgb_img)


% Compress/rescale the irt image
[rows, columns, numberOfColorChannels] = size(irt_img_bmp);
newWidth = round(0.87 * columns)
stretchedImage = imresize(irt_img_bmp, [rows newWidth]);
irt_img_bmp = stretchedImage;


% Uncomment this if want to select the movingPoints and fixedPoints
% leg_fixed = rgb_img;
% leg_moving = irt_img; 
% cpselect(leg_moving, leg_fixed);


% % treshold the irt img to highlight the high gradient region
% [bw, rgbImage] = irtThresholder(irt_img_bmp);
% irt_img_bmp = rgbImage;



% % display both images
% figure
% imshowpair(rgb_img, irt_img_bmp, 'montage');


% Define some matching control points on the fixed image (RGB img) and
% moving image (IRT img) ontain from 
fixedPoints = [6.031250000000001e+02 2.461249999999999e+02;8.903750000000002e+02 1.103749999999999e+02;6.128750000000002e+02 6.181250000000000e+02;7.688750000000001e+02 6.361250000000000e+02];
movingPoints = [1.531250000000000e+02 2.078749999999999e+02;4.328750000000001e+02 22.624999999999886;2.348750000000001e+02 7.343750000000000e+02;4.021250000000001e+02 7.501250000000000e+02];


% Create a geometric transformation that can be used to align the two images
tform = fitgeotrans(movingPoints,fixedPoints,'similarity')


% Use the tform estimate to resample the rotated image to register it with the fixed image. 
% Overlap the irt image on the rgb image
irt_registered = imwarp(irt_img,tform,'OutputView',imref2d(size(rgb_img)));
figure
imshowpair(rgb_img,irt_registered,'blend')

