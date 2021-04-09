function [data] = handling_transient(file,varargin)
%HANDLING_TRANSIENT takes in .csv transient data and processes it into a
%form consistent with objective function output.
%   HANDLING_TRANSIENT takes the name of a .csv file as input. The .csv
%   file should be formatted as a table with headings and column-wise data.
%   The case-insensitive headings of "time" "vgs" 'vds" and "id" are
%   extracted and returned along with a heading cell array.

varg = varargin{:};

if ~isempty(varg)
    s = varg{1};
else
    s = 1;
end

if numel(varg)>1
    targ_heads = varg{2};
else
    targ_heads = {'time', 'vds', 'vgs', 'id'};
end

f1=strcat("Fitting_Data\Transient_Data\",file);
raw_table = readtable(f1);

raw_heads = raw_table.Properties.VariableNames;

[~, ia, ib] = intersect(lower(raw_heads), targ_heads);

[~, is] = sort(ib);

ia = ia(is);
ib = ib(is);

data = targ_heads(ib);

for k = 1:size(ib,1)
    cur_data = raw_table.(raw_heads{ia(k)});
    if ~strcmpi('time', raw_heads{ia(k)})
        cur_data = smooth(cur_data,s);
    end
    if iscolumn(cur_data)
        cur_data = cur_data';
    end
    data{2,k} = cur_data;
end

% headings = targ_heads(ib);
% 
% data = [];
% for k = 1:size(ib,1)
%     cur_data = raw_table.(raw_heads{ia(k)});
%     if ~strcmpi('time', raw_heads{ia(k)})
%         cur_data = smooth(cur_data,s);
%     end
%     data = [data,cur_data];
% end

end
