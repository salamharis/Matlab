close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img = 'IRT_31_r';
rgb_img = 'RGB_31_r';

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

% imshow(irt_img)

% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));


[bw,rgb_threshold] = legThreshold3(rgb_img);

bw = bwareaopen(bw,100);
bw = imfill(bw,'holes');
% Find edges using the Canny operator with hysteresis thresholds of 0.1
% and 0.2 with smoothing parameter sigma set to 1.
edgeim_rgb = edge(bw,'canny', [0.1 0.2], 1);
edgeim_irt = edge(irt_img,'canny', [0.1 0.2], 1);

figure(1), imshow(edgeim_rgb);
title('RGB Image Edges');
figure(2), imshow(edgeim_irt);
title('IRT Image Edges');


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
tform_mean = [0.878866245374861 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;3.605920167976121e+02 52.808982154224054 1];
tform_mean = affine2d(tform_mean);

irt_registered = imwarp(edgeim_irt,tform_mean,'OutputView',imref2d(size(rgb_img)));
% edgeim_irt = edge(irt_registered,'canny', [0.1 0.2], 1);
% img_title = ['Image ',num2str(p),' position ',leg];
figure(3), imshow(irt_registered);
title('IRT Image Edges Registered');
figure(4), imshowpair(rgb_img,irt_registered,'blend')
title('IRT and RGB Image Edges Registered');
% title(img_title);

width_rgb = numel(bw(1,:));
height_rgb = numel(bw(:,1));
width_irt = numel(irt_registered(1,:));
height_irt = numel(irt_registered(:,1));


% Detect the line coordinates on the image
y_line_rgb{1,1} = height_rgb/3;
y_line_rgb{1,2} = height_rgb/2;
y_line_rgb{1,3} = height_rgb*2/3;
horizontal_line_rgb = zeros(1,width_rgb);
index = 1;
for i=1 : width_rgb
    if edgeim_rgb(y_line_rgb{1,1},i) == 1
    horizontal_line_rgb(index) = i;
    index = index + 1;
    end
end


y_line_irt{1,1} = height_irt/3;
y_line_irt{1,2} = height_irt/2;
y_line_irt{1,3} = height_irt*2/3;
horizontal_line_irt = zeros(1,width_irt);
index = 1;
for i=1 : width_irt
    if irt_registered(y_line_irt{1,1},i) == 1
    horizontal_line_irt(index) = i;
    index = index + 1;
    end
end

% to ensure there no near edge / to reduce false edge
num_edge_point_rgb = nnz(horizontal_line_rgb);
for i=2 : num_edge_point_rgb 
    k = abs(horizontal_line_rgb(i)-horizontal_line_rgb(i-1)) < 10;
    if k == 1
        horizontal_line_rgb(i-1) = abs(floor((horizontal_line_rgb(i-1)+horizontal_line_rgb(i))/2));
        % shift the array
        for ii=i : num_edge_point_rgb
            horizontal_line_rgb(ii) = horizontal_line_rgb(ii+1);        
        end
        % make sure the last array not extra value after shifted 
        horizontal_line_rgb(num_edge_point_rgb) = 0;
        num_edge_point_rgb = nnz(horizontal_line_rgb); 
    end
end

num_edge_point_irt = nnz(horizontal_line_irt);
for i=2 : num_edge_point_irt 
    k = abs(horizontal_line_irt(i)-horizontal_line_irt(i-1)) < 10;
    if k == 1
        horizontal_line_irt(i-1) = abs(floor((horizontal_line_irt(i-1)+horizontal_line_irt(i))/2));
        % shift the array
        for ii=i : num_edge_point_irt
            horizontal_line_irt(ii) = horizontal_line_irt(ii+1);        
        end
        % make sure the last array not extra value after shifted 
        horizontal_line_irt(num_edge_point_irt) = 0;
        num_edge_point_irt = nnz(horizontal_line_irt); 
    end
end


% move the irt image to the rgb image location
edge_diff = horizontal_line_rgb(1,1) - horizontal_line_irt(1,1);
new_translation_x = tform_mean.T(3,1) + edge_diff;
tform_mean = [0.878866245374861 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;new_translation_x 52.808982154224054 1];
tform_mean = affine2d(tform_mean);

irt_registered = imwarp(edgeim_irt,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(5), imshowpair(rgb_img,irt_registered,'blend')
title('After Shift the IRT Image')
% irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
% figure(6), imshowpair(rgb_img,irt_registered,'blend')
% title('Final Output')


% calculate the distance between the leg edges for measure leg diameter
leg_diameter_rgb = horizontal_line_rgb(1,2) - horizontal_line_rgb(1,1);
leg_diameter_irt = horizontal_line_irt(1,2) - horizontal_line_irt(1,1);
diff_diameter = leg_diameter_rgb - leg_diameter_irt;
if diff_diameter>20  %rgb wider than irt
    % scale the x axis irt wider
    new_scale_x = tform_mean.T(1,1)*(leg_diameter_irt/leg_diameter_rgb);
elseif diff_diameter<20 %irt wider than irt
    % scale the x axis irt narrow
    new_scale_x = tform_mean.T(1,1)*(leg_diameter_rgb/leg_diameter_irt);
else
end

tform_mean = [new_scale_x 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;new_translation_x 52.808982154224054 1];
tform_mean = affine2d(tform_mean);
irt_registered = imwarp(edgeim_irt,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(6), imshowpair(rgb_img,irt_registered,'blend')
title('After Scale the IRT Image')

