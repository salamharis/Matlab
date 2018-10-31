path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img = 'IRT_19_l';
rgb_img = 'RGB_19_l';

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
imshow(irt_img)

% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));


% Compress/rescale the irt image
[rows, columns, numberOfColorChannels] = size(irt_img);
newWidth = round(0.7734 * columns);
stretchedImage = imresize(irt_img, [rows newWidth]);
irt_img = stretchedImage;


% Uncomment this if want to select the movingPoints and fixedPoints
leg_fixed = rgb_img;
leg_moving = irt_img; 
cpselect(leg_moving, leg_fixed);



% % treshold the irt img to highlight the high gradient region
% [bw, rgbImage] = irtThresholder(irt_img_bmp);
% irt_img_bmp = rgbImage;



% % display both images
% figure
% imshowpair(rgb_img, irt_img_bmp, 'montage');


% Define some matching control points on the fixed image (RGB img) and
% moving image (IRT img) ontain from 
% fixedPoints = [7.606250000000001e+02 1.201249999999999e+02;9.616250000000002e+02 48.874999999999886;7.471250000000001e+02 4.411250000000000e+02;8.648750000000002e+02 5.093750000000000e+02;9.773750000000002e+02 2.498749999999999e+02];
% movingPoints = [2.172500000000000e+02 57.749999999999886;3.717500000000000e+02 16.249999999999886;2.393421052631580e+02 3.577293233082707e+02;3.222500000000000e+02 4.342500000000000e+02;3.892500000000000e+02 2.042499999999999e+02];


% % Create a geometric transformation that can be used to align the two images
% tform = fitgeotrans(movingPoints,fixedPoints,'affine')
% 
% 
% % Use the tform estimate to resample the rotated image to register it with the fixed image. 
% % Overlap the irt image on the rgb image
% irt_registered = imwarp(irt_img,tform,'OutputView',imref2d(size(rgb_img)));
% figure
% imshowpair(rgb_img,irt_registered,'montage')