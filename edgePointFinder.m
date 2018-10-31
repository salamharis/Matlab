function [edge_point,num_edge_point] = edgePointFinder(edge_img,y_line,width_img)

% edge_point = zeros(1,width_img);
% index = 1;
% for i=1 : width_img
%     if edge_img(y_line,i) == 1
%     edge_point(index) = i;
%     index = index + 1;
%     end
% end

edge_point = find(edge_img(y_line,:)==1);

% to ensure there no near edge / to reduce false edge
num_edge_point = nnz(edge_point);
% for i=2 : num_edge_point 
%     k = abs(edge_point(i)-edge_point(i-1)) < 10;
%     if k == 1
%         edge_point(i-1) = abs(floor((edge_point(i-1)+edge_point(i))/2));
%         % shift the array
%         for ii=i : num_edge_point
%             edge_point(ii) = edge_point(ii+1);        
%         end
%         % make sure the last array not extra value after shifted 
%         edge_point(num_edge_point) = 0;
%         num_edge_point = nnz(edge_point); 
%     end
% end

end

