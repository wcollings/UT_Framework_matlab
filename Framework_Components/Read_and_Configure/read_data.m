function [ m_data ] = read_data( csv_file, varargin )
%READDATA reads data from a .csv. If time data is present, it is shifted to
%begin at 0. Other data is smoothed with a span of s if provided.  
%   Input .csv file is expected in the form of a column table with
%   headings. All headings are made lower case for consistency. Output 
%   m_data is a cell array in which the first row of cells contains column 
%   headings and the second row of cells contains column vectors of data 
%   corresponding to the heading. 

%%%
% If there is a variable argument given, the first value provided is a
% smoothing span for a moving average of the data. 
s = 0; 
if nargin > 1
    s = varargin{1};
end

%%%
% Extract data from the .csv file and prepare for formatting. Create a
% table from the .csv data. The number of headings corresponds to the
% number of variables given. The output value m_data is initialized with
% this knowledge for speed. 
data = readtable(csv_file);
headings = data.Properties.VariableNames;
n_vars = numel(headings);
m_data = cell(2, n_vars);

%%%
% Properly format the extracted data and insert it into the already
% initialized cell array. If the case-insensitive heading 'time' is
% present, the corresponding data will be shifted to begin at zero. If a
% smoothing span was provided, take a moving average of the data with this
% span. Format the output data with the first row of the cell array
% containing the heading names (lower case for consistency) and the second
% row containing row vectors. 
for i = 1:n_vars
    cur_data = data{:,i};
    if strcmpi(headings{1,i}, 'time')
        t_shift = cur_data(1);
        cur_data = cur_data - t_shift;
    else
        if s
            cur_data = smooth(cur_data, s, 'moving');
        end
    end
    if iscolumn(cur_data)
        cur_data = cur_data'; 
    end
    m_data{2,i} = cur_data;
end

m_data(1,:) = lower(headings);
end

