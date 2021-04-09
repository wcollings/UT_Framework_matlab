function [approximated_data] = RW(data,epsilon)
%RW takes a dataset and epsilon value and produces the
%Reumann-Witkam approximation of the data.
%   data is an mxn dataset with dimensionality m and number of points n.

if size(data,1) == 2
    data = [data;zeros(1,size(data,2))];
end

% Identify first and last points as keys (to be kept)
keys = data(:,1);
p_end = data(:,end);

k = 2;
while k < size(data,2)
    p1 = keys(end);
    p2 = data(:,k);
    remaining_data = data(:,k:end);
    a = p1 - p2;
    b = remaining_data - p2;
    dists = zeros(1,size(data,2));
    for r = 1:size(b,2)
        dists(r+1) = norm(cross(a,b(:,r)))/norm(a);
    end
    miss_idx = find(dists>epsilon,1);
    new_key = remaining_data(:,miss_idx-1);
    keys = [keys, new_key];
    k = k + miss_idx;
end

approximated_data = [keys, p_end];

end

