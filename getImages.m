path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img = 'IRT_01_l';
rgb_img = 'RGB_01_l';

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

% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));

% display both images
figure
imshowpair(rgb_img, irt_img, 'montage');

%Define some matching control points on the fixed image (RGB img) and moving image (IRT img).
fixedPoints = [1121 264; 937 596];
movingPoints = [758 264; 584 596];

%Create a geometric transformation that can be used to align the two images
tform = fitgeotrans(movingPoints,fixedPoints,'NonreflectiveSimilarity')

%Use the tform estimate to resample the rotated image to register it with the fixed image. 
%The regions of color (green and magenta) in the false color overlay image indicate error in the registration. 
%This error comes from a lack of precise correspondence in the control points.
irt_registered = imwarp(irt_img,tform,'OutputView',imref2d(size(rgb_img)));
figure
imshowpair(rgb_img,irt_registered)
