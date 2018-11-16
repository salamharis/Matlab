function [score1,score2] = bestLeg(leg_value1,leg_value2,x_midpoint_irt)
score1 = 0; %first leg score
score2 = 0; %second leg score
if leg_value2~=0 %to ensure there are second leg
    dist_from_mid1 = abs(x_midpoint_irt - (leg_value1(2)+leg_value1(1))/2);
    dist_from_mid2 = abs(x_midpoint_irt - (leg_value2(2)+leg_value2(1))/2);
    if dist_from_mid1<dist_from_mid2
        score1=score1+1;
    else
        score2=score2+1;
    end
else
    score1=score1+1;
end
end

