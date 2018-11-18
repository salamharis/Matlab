close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path = '\\134.130.86.237\projekt\vulnusMON\201802_Bochum\aufnahmen';
path = 'C:\Users\MeoW\Desktop\Image Registration\Picture-20181016T194623Z-001\Picture';
addpath(path);
img_num = '26';
leg_position = 'l';

irt_img_name = ['IRT_',img_num,'_',leg_position];
rgb_img_name = ['RGB_',img_num,'_',leg_position];

% get 2D matrix (irt_img) containing temperature values for each pixel in °C
irt_img = dir(strcat(path,'\',irt_img_name,'*.asc'));
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
rgb_img = imread(strcat(path,'\',rgb_img_name,'.jpg'));
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

%determine the group of transformation matrix
img_num = str2double(img_num);
img_group = 1*(img_num>=3&&img_num<=21)+2*(img_num>=25&&img_num<=27)+3*(img_num>=28&&img_num<=32)...
            +4*(img_num>=33&&img_num<=38);
if img_group==1
 %4-21
 tform_mean = [0.805599525356956 0.020598075940177 0;-0.104893577450360 0.781929514272011 0;4.302242013663251e+02 48.320819733806360 1];
elseif img_group==2
 %25-27
 tform_mean = [0.637487839738744 -0.004434008430408 0;-0.040811608158753 0.551993295992885 0;4.416109116294455e+02 1.549737714257566e+02 1];
elseif img_group==3
 %28-32
tform_mean = [0.878866245374861 0.009763143750244 0;-0.065736746411131 0.788104250865920 0;3.605920167976121e+02 52.808982154224054 1];
elseif img_group==4
%33-38
 tform_mean =[0.822833383351964 8.876458614865464e-04 0;-0.023691893394736 0.829565435596974 0;2.270445581698194e+02 56.856631985064155 1];
else
    %if no group, just use the most used transformation matrix
    tform_mean = [0.805599525356956 0.020598075940177 0;-0.104893577450360 0.781929514272011 0;4.302242013663251e+02 48.320819733806360 1];
end
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
if numel(edgepoint_upper_rgb)>1 && numel(edgepoint_upper_irt)>1 %to make sure the are sufficient edge
    first_edge_diff = abs(edgepoint_upper_rgb(1,1) - edgepoint_upper_irt(1,1));
    second_edge_diff = abs(edgepoint_upper_rgb(1,2) - edgepoint_upper_irt(1,2));
    %if the differences more than 10 pixels
    pixel_diff = 10;
    additional_adjust = first_edge_diff>pixel_diff & second_edge_diff>pixel_diff;
else
    additional_adjust = false;
end
additional_adjust = true; %///testing///, comment when not use
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

min_leg_gradient = 0.5;
g = irt_img.*(irt_img>min_leg_gradient);
figure(5),imshow(g);
leg = cell(4,3);
leg = {'Line' 'first_leg_range' 'second_leg_range';'Upper' 0 0;'Mid' 0 0;'Bottom' 0 0};

[leg{2,2},leg{2,3}]=irtEdge(irt_img,y_line_upper_irt,min_leg_gradient);
[leg{3,2},leg{3,3}]=irtEdge(irt_img,y_line_mid_irt,min_leg_gradient);
[leg{4,2},leg{4,3}]=irtEdge(irt_img,y_line_bottom_irt,min_leg_gradient);


% Choosing the leg, the most centered and highest irt intensity
x_midpoint_irt = irt_width/2;
score_leg1 = 0; %first leg score
score_leg2 = 0; %second leg score
[score1,score2] = bestLeg(leg{2,2},leg{2,3},x_midpoint_irt);
score_leg1=score_leg1+score1;
score_leg2=score_leg2+score2;
[score1,score2] = bestLeg(leg{3,2},leg{3,3},x_midpoint_irt);
score_leg1=score_leg1+score1;
score_leg2=score_leg2+score2;
[score1,score2] = bestLeg(leg{4,2},leg{4,3},x_midpoint_irt);
score_leg1=score_leg1+score1;
score_leg2=score_leg2+score2;


% Find if the leg have reasonable diameter
if score_leg1>score_leg2
    %find the shortest diameter/smallest x-axis point
    compare_point = 1*(leg{2,2}(2)<leg{3,2}(2)&leg{2,2}(2)<leg{4,2}(2))+...
        2*(leg{3,2}(2)<leg{2,2}(2)&leg{3,2}(2)<leg{4,2}(2))+...
        3*(leg{4,2}(2)<leg{2,2}(2)&leg{4,2}(2)<leg{3,2}(2));
    x_small_point = (compare_point==1)*leg{2,2}(2)+(compare_point==2)*leg{3,2}(2)+(compare_point==3)*leg{4,2}(2);
    y_small_point = (compare_point==1)*y_line_upper_irt+(compare_point==2)*y_line_mid_irt+(compare_point==3)*y_line_bottom_irt;
    
    
    size_scan = 6; %size of adjacent pt left and right
    upward_dist = y_small_point-1;
    right_small_pt = ones(upward_dist,size_scan);
    left_small_pt = ones(upward_dist,size_scan);
    im_newline = ones(irt_height, irt_width);
    for f=1:upward_dist
        
        % find the point that in between two gradient
        best_edge = false;
        while ~best_edge
            
            %find the pixel value in horizontal and adjacent value
            right_small_pt(f,:) = irt_img(y_small_point,x_small_point+1:x_small_point+size_scan);
            left_small_pt(f,:) = irt_img(y_small_point,x_small_point-size_scan:x_small_point-1);
            right_small_diff = [diff(right_small_pt(f,:)),0];
            left_small_diff = [diff(left_small_pt(f,:)),0];
            mean_right = abs(mean(right_small_diff));
            mean_left = abs(mean(left_small_diff));
            mean_thres = 0.01;
            
            if mean_left>0&&mean_right>0 %when there are heat gradient at both side
                
                if mean_right>mean_thres&&mean_left>mean_thres  %the pt is between two leg
                    %do nothing, the point is at edge
                    best_edge = true;
                elseif mean_right>mean_thres&&mean_left<mean_thres  %the pt at the leg area
                    %shift right
                    x_small_point=x_small_point+1;
                elseif mean_right<mean_thres&&mean_left>mean_thres  %the pt at the outer leg area
                    %shift left
                    x_small_point=x_small_point-1;
                elseif mean_right<mean_thres&&mean_left<mean_thres  %the pt is far away from the leg
                    best_edge = true;
                end
            else
                
                if any(~left_small_pt(f,:))&&all(~right_small_pt(f,:)) %if there are zero value on left
                    zero_pt = find(~(left_small_pt(f,:)));
                    x_small_point=x_small_point-(size_scan+1-zero_pt(1,1)); %shift to left
                    
                elseif all(left_small_pt(f,:))&&any(right_small_pt(f,:))
                    one_pt = find(right_small_pt(f,:));
                    x_small_point=x_small_point+one_pt(1,max(one_pt)); %shift to right
                    
                else
                    %do nothing, the point is at edge
                    best_edge = true;
                end
            end
        end   
        
        %shift the point upward about 1
        y_small_point = y_small_point-1;
        
        %draw the new edge
        im_newline(y_small_point,x_small_point)=0;
    end
    figure(6), imshow(im_newline);
    figure(7), imshowpair(irt_img,im_newline,'blend');
    
end
irt_img_new = irt_img.*im_newline;
[leg{2,2},leg{2,3}]=irtEdge(irt_img_new,y_line_upper_irt,min_leg_gradient);
[leg{3,2},leg{3,3}]=irtEdge(irt_img_new,y_line_mid_irt,min_leg_gradient);
[leg{4,2},leg{4,3}]=irtEdge(irt_img_new,y_line_bottom_irt,min_leg_gradient);

max_intensity = 0.7;

% leg_diameter_rgb = edgepoint_upper_rgb(1,2) - edgepoint_upper_rgb(1,1);
% leg_diameter_irt = edgepoint_upper_irt(1,2) - edgepoint_upper_irt(1,1);
leg_diameter_rgb = edgepoint_upper_rgb(1,2) - edgepoint_upper_rgb(1,1);
leg_diameter_irt = leg{2,2}(2) - leg{2,2}(1);
diff_diameter = leg_diameter_rgb - leg_diameter_irt;
rgb_wider = diff_diameter>20;
irt_wider = diff_diameter<20;

a = tform_mean.T(1,1);  %read the previous transform matrix
b = tform_mean.T(1,2);
c = tform_mean.T(2,1);
d = tform_mean.T(2,2);
e = tform_mean.T(3,1);
f = tform_mean.T(3,2);

a = rgb_wider*a*(leg_diameter_irt/leg_diameter_rgb)+... % scale the x axis irt wider
    irt_wider*a*(leg_diameter_rgb/leg_diameter_irt);    % scale the x axis irt narrow


tform_mean = [a b 0;c d 0;e f 1];
tform_mean = affine2d(tform_mean);
irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(8), imshowpair(rgb_img,irt_registered,'blend')
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
edgepoint_mid_irt = find(edgeim_irt(y_line_upper_irt,:)==1);
edge_diff = edgepoint_mid_rgb(1,1) - edgepoint_mid_irt(1,1);
e = e + edge_diff;
tform_mean = [a b 0;c d 0;e f 1];
tform_mean = affine2d(tform_mean);
irt_registered = imwarp(irt_img,tform_mean,'OutputView',imref2d(size(rgb_img)));
figure(9), imshowpair(rgb_img,irt_registered,'blend')
title('After Shift the IRT Image')

end


