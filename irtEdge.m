function [leg_value1,leg_value2]=irtEdge(irt_img,y_line,min_leg_gradient)
line_irt = irt_img(y_line,:); %find the horinzontal value
leg_area = find(line_irt.*(line_irt>min_leg_gradient)); %to remove other than leg information
area_diff = [diff(leg_area),0];
leg_value1(1) = leg_area(1);
temp = leg_area(area_diff>1);  %to find the first leg
if any(temp)&&temp(1)~=max(leg_area) % Determine if there are second leg
    temp2 = temp(1)==leg_area;
    temp2 = [temp2(end),temp2(1:end-1)]; %to find the second leg point
    temp2 = leg_area(temp2);
    leg_value1(2) = temp(1);
    leg_value2(1) = temp2;
    leg_value2(2) = max(leg_area);  %end point of second leg
else
    leg_value1(2) = max(leg_area);  %first leg end point
    leg_value2 = 0; %there no second leg, or two leg overlaping
end

end
