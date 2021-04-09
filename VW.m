function [approximated_data] = VW(data,epsilon)
%VW takes a dataset and epsilon value and produces the
%Visvalingam & Whyatt approximation of the data.
%   data is an mxn dataset with dimensionality m and number of points n.

if size(data,1) == 2
    data = [data;zeros(1,size(data,2))];
end

approximated_data = data; 
areas = NaN(1,size(data,2)-2);

while any(areas<epsilon|isnan(areas))
%     if ~isnan(areas)
%         small_areas = find(areas<epsilon);
%         areas(small_areas)=[];
%         approximated_data(:,small_areas+1)=[];
%     end
    [~, min_idx] = min(areas);
    approximated_data(:,min_idx+1) = [];
    areas(min_idx) = [];
    for k = 1:size(areas,2)
        p1 = approximated_data(:,k);
        p2 = approximated_data(:,k+1);
        p3 = approximated_data(:,k+2);
        v1 = p2-p1; 
        v2 = p3-p1;
        areas(k) = norm(cross(v1,v2));
    end
    
end



end

