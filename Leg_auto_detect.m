close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
addpath(path);

irt_img = 'IRT_30_r';
rgb_img = 'RGB_30_r';

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

% irt_img = irt_img{1,1}(16:size(irt_img{1,1}),1);
% irt_img = str2double(strrep(irt_img,',','.'));
% irt_img = reshape(irt_img,[irt_width,irt_height]);
% irt_img = irt_img';

irt_img = (irt_img - min(irt_img(:)))/(max(irt_img(:)) - min(irt_img(:)));
% imshow(irt_img)

% get 3D matrix (rgb_img) containing R, G and B color values for each pixel
rgb_img = imread(strcat(path,'\',rgb_img,'.jpg'));
rgb_img = imresize(rgb_img,min(irt_width,irt_height)/min(size(rgb_img(:,:,1))));
min_leg_gradient = 0.5;
irt_img = irt_img.*(irt_img>min_leg_gradient); %remove unwanted background value

[bw,rgb_threshold] = legThreshold3(rgb_img);

bw = bwareaopen(bw,100);
bw = imfill(bw,'holes');
% Find edges using the Canny operator with hysteresis thresholds of 0.1
% and 0.2 with smoothing parameter sigma set to 1.
edgeim_rgb = edge(bw,'canny', [0.1 0.2], 10);
edgeim_irt = edge(irt_img,'canny', [0.1 0.2], 10);

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
%28-33
tform_mean = [0.878866245374861 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;3.605920167976121e+02 52.808982154224054 1];
%33-38
% tform_mean =[0.822833383351964 8.876458614865464e-04 0;-0.023691893394736 0.829565435596974 0;2.270445581698194e+02 56.856631985064155 1];
tform_mean = affine2d(tform_mean);

irt_registered = imwarp(edgeim_irt,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(3), imshowpair(rgb_img,irt_registered,'blend')
title('IRT and RGB Image Edges Registered');
irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(4), imshowpair(rgb_img,irt_registered,'blend')
title('IRT and RGB Image Registered');


width_rgb = numel(bw(1,:));
height_rgb = numel(bw(:,1));
width_irt = numel(irt_registered(1,:));
height_irt = numel(irt_registered(:,1));


% Detect the line coordinates on the image on 
% three different horizontal line
y_line_upper_rgb = height_rgb/3;
y_line_mid_rgb = height_rgb/2;
y_line_bottom_rgb = height_rgb*2/3;

y_line_upper_irt = height_irt/3;
y_line_mid_irt = height_irt/2;
y_line_bottom_irt = height_irt*2/3;


% Find the edge point
edgeim_irt = edge(irt_registered,'canny', [0.1 0.2], 10);

edgepoint_upper_rgb = find(edgeim_rgb(y_line_upper_rgb,:)==1);
edgepoint_mid_rgb = find(edgeim_rgb(y_line_mid_rgb,:)==1);
edgepoint_bottom_rgb = find(edgeim_rgb(y_line_bottom_rgb,:)==1);

edgepoint_upper_irt = find(edgeim_irt(y_line_upper_irt,:)==1);
edgepoint_mid_irt = find(edgeim_irt(y_line_mid_irt,:)==1);
edgepoint_bottom_irt = find(edgeim_irt(y_line_bottom_irt,:)==1);

% Determine does the image required additional adjustment
first_edge_diff = abs(edgepoint_upper_rgb(1,1) - edgepoint_upper_irt(1,1));
second_edge_diff = abs(edgepoint_upper_rgb(1,2) - edgepoint_upper_irt(1,2));
%if the differences more than 10 pixels
pixel_diff = 10;
additional_adjust = first_edge_diff>pixel_diff & second_edge_diff>pixel_diff; 

if additional_adjust == true
%% Simplify the transformation function
%     [a   b   0;
%      c   d   0;
%      e   f   1]
a = tform_mean.T(1,1);
b = tform_mean.T(1,2);
c = tform_mean.T(2,1);
d = tform_mean.T(2,2);
e = tform_mean.T(3,1);
f = tform_mean.T(3,2);


%% Scale the irt image
% calculate the distance between the leg edges for measure leg diameter
% edgeim_irt = edge(irt_registered,'canny', [0.1 0.2], 10);
% edgepoint_upper_irt = find(edgeim_irt(y_line_upper_irt,:)==1);
% edgepoint_mid_irt = find(edgeim_irt(y_line_mid_irt,:)==1);
% edgepoint_bottom_irt = find(edgeim_irt(y_line_bottom_irt,:)==1);


% find the gradient of the leg
upperline_irt = irt_img(y_line_upper_irt,:);
min_leg_gradient = 0.5;
g = irt_img.*(irt_img>min_leg_gradient);
figure(7),imshow(g);
leg_area = find(upperline_irt.*(upperline_irt>min_leg_gradient));
area_diff = [diff(leg_area),0];
leg = cell(4,3);
leg = {'Line' 'first_leg_range' 'second_leg_range';'Upper' 0 0;'Mid' 0 0;'Bottom' 0 0};
leg{2,2}(1) = leg_area(1);
temp = leg_area(area_diff>2);  %to find the first leg
if any(temp)&&temp(1)~=max(leg_area) % Determine if there are second leg
    temp2 = temp(1)==leg_area;
    temp2 = [temp2(end),temp2(1:end-1)]; %to find the second leg point
    temp2 = leg_area(temp2);
    leg{2,2}(2) = temp(1);
    leg{2,3}(1) = temp2;
    leg{2,3}(2) = max(leg_area);  %end point of second leg
else
    leg{2,2}(2) = max(leg_area);  %first leg end point
end

midline_irt = irt_img(y_line_mid_irt,:);
leg_area = find(midline_irt.*(midline_irt>min_leg_gradient));
area_diff = [diff(leg_area),0];
leg{3,2}(1) = leg_area(1);
temp = leg_area(area_diff>2);
if any(temp)&&temp(1)~=max(leg_area)
    temp2 = temp(1)==leg_area;
    temp2 = [temp2(end),temp2(1:end-1)]; %to find the second leg point
    temp2 = leg_area(temp2);
    leg{3,2}(2) = temp(1);
    leg{3,3}(1) = temp2;
    leg{3,3}(2) = max(leg_area);  %end point of second leg
else
    leg{3,2}(2) = max(leg_area);
end

bottomline_irt = irt_img(y_line_bottom_irt,:);
leg_area = find(bottomline_irt.*(bottomline_irt>min_leg_gradient));
area_diff = [diff(leg_area),0];
leg{4,2}(1) = leg_area(1);
temp = leg_area(area_diff>2);
if any(temp)&&temp(1)~=max(leg_area)
    temp2 = temp(1)==leg_area;
    temp2 = [temp2(end),temp2(1:end-1)]; %to find the second leg point
    temp2 = leg_area(temp2);
    leg{4,2}(2) = temp(1);
    leg{4,3}(1) = temp2;
    leg{4,3}(2) = max(leg_area);  %end point of second leg
else
    leg{4,2}(2) = max(leg_area);
end

% Choosing the leg, the most centered and highest irt intensity
x_midpoint_irt = irt_width/2;
score1 = 0; %first leg score
score2 = 0; %second leg score
if leg{2,3}~=0 %to ensure there are second leg
    dist_from_mid1 = abs(x_midpoint_irt - (leg{2,2}(2)+leg{2,2}(1))/2);
    dist_from_mid2 = abs(x_midpoint_irt - (leg{2,3}(2)+leg{2,3}(1))/2);
    if dist_from_mid1<dist_from_mid2
        score1=score1+1;
    else
        score2=score2+1;
    end
else
    score1=score1+1;
end
if leg{3,3}~=0
    dist_from_mid1 = abs(x_midpoint_irt - (leg{3,2}(2)+leg{3,2}(1))/2);
    dist_from_mid2 = abs(x_midpoint_irt - (leg{3,3}(2)+leg{3,3}(1))/2);
    if dist_from_mid1<dist_from_mid2
        score1=score1+1;
    else
        score2=score2+1;
    end
else
    score1=score1+1;
end
if leg{4,3}~=0
    dist_from_mid1 = abs(x_midpoint_irt - (leg{4,2}(2)+leg{4,2}(1))/2);
    dist_from_mid2 = abs(x_midpoint_irt - (leg{4,3}(2)+leg{4,3}(1))/2);
    if dist_from_mid1<dist_from_mid2
        score1=score1+1;
    else
        score2=score2+1;
    end
else
    score1=score1+1;
end

% Find if the leg have reasonable diameter
if score1>score2
    %find the shortest diameter/smallest x-axis point
    compare_point = 1*(leg{2,2}(2)<leg{3,2}(2)&leg{2,2}(2)<leg{4,2}(2))+...
                  2*(leg{3,2}(2)<leg{2,2}(2)&leg{3,2}(2)<leg{4,2}(2))+...
                  3*(leg{4,2}(2)<leg{2,2}(2)&leg{4,2}(2)<leg{3,2}(2));
    x_small_point = (compare_point==1)*leg{2,2}(2)+(compare_point==2)*leg{3,2}(2)+(compare_point==3)*leg{4,2}(2);
    
    %find the pixel value in horizontal and adjacent value
    y_small_point = (compare_point==1)*y_line_upper_irt+(compare_point==2)*y_line_mid_irt+(compare_point==3)*y_line_bottom_irt;
    y_small_point = y_small_point-1;
    right_small_pt = [irt_img(y_small_point,x_small_point+1),irt_img(y_small_point,x_small_point+2),irt_img(y_small_point,x_small_point+3)]
    left_small_pt = [irt_img(y_small_point,x_small_point-1),irt_img(y_small_point,x_small_point-2),irt_img(y_small_point,x_small_point-3)]
end



max_intensity = 0.7;

leg_diameter_rgb = edgepoint_upper_rgb(1,2) - edgepoint_upper_rgb(1,1);
leg_diameter_irt = edgepoint_upper_irt(1,2) - edgepoint_upper_irt(1,1);
diff_diameter = leg_diameter_rgb - leg_diameter_irt;
rgb_wider = diff_diameter>20;
irt_wider = diff_diameter<20;
a = rgb_wider*a*(leg_diameter_irt/leg_diameter_rgb)+... % scale the x axis irt wider
    irt_wider*a*(leg_diameter_rgb/leg_diameter_irt);    % scale the x axis irt narrow


tform_mean = [a b 0;c d 0;e f 1];
tform_mean = affine2d(tform_mean);
irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(5), imshowpair(rgb_img,irt_registered,'blend')
title('After Scale the IRT Image')

% %% Rotation transformation function
% %     [cos(q)    sin(q)   0;
% %      -sin(q)   cos(q)   0;
% %      0         0        1]
% 
% % measure the angle need to rotate
% edgeim_irt = edge(irt_registered,'canny', 0.2, 1);
% edgepoint_upper_irt = find(edgeim_irt(y_line_upper_irt,:)==1);
% edgepoint_bottom_irt = find(edgeim_irt(y_line_bottom_irt,:)==1);
% 
% length_theta_x_irt = edgepoint_upper_irt(1,2) - edgepoint_bottom_irt(1,2);
% length_theta_y_irt = y_line_bottom_irt - y_line_upper_irt;
% length_theta_x_rgb = edge_point_upper_rgb(1,1) - edge_point_bottom_rgb(1,1);
% length_theta_y_rgb = y_line_bottom_rgb - y_line_upper_rgb;
% 
% theta_irt = atand(length_theta_x_irt/length_theta_y_irt);
% theta_rgb = atand(length_theta_x_rgb/length_theta_y_rgb);
% theta = theta_rgb - theta_irt;
% min_angle_diff = 5;
% if abs(theta)>min_angle_diff
% rotate_matrix = [cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0;0 0 1];
% tform_mean = rotate_matrix*tform_mean.T;
% tform_mean = affine2d(tform_mean);
% irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
% figure(7), imshowpair(rgb_img,irt_registered,'blend')
% title('After rotate the IRT Image')
% end

%% Shift the irt image to the rgb image location
a = tform_mean.T(1,1);
b = tform_mean.T(1,2);
c = tform_mean.T(2,1);
d = tform_mean.T(2,2);
e = tform_mean.T(3,1);
f = tform_mean.T(3,2);

edgeim_irt = edge(irt_registered,'canny', [0.1 0.2], 10);
edgepoint_upper_irt = find(edgeim_irt(y_line_upper_irt,:)==1);
edge_diff = edgepoint_upper_rgb(1,1) - edgepoint_upper_irt(1,1);
e = e + edge_diff;
tform_mean = [a b 0;c d 0;e f 1];
tform_mean = affine2d(tform_mean);
irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(6), imshowpair(rgb_img,irt_registered,'blend')
title('After Shift the IRT Image')

end



