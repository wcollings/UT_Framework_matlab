function [data] = handling_fwdcurves(data_file, varargin)
%%%HANDLING_FWDCURVES takes a data .csv file in which data is organized in
%columns with headings along with any handling arguments provided and
%returns processed data from the data file.
%   Data is output as a 1x3 cell array with three entries: the first cell
%   contains a vector of gate voltages at which forward curves have been
%   taken, the second contains a vector of drain voltages, and the third
%   contains a matrix of drain currents where each column corresponds to an
%   individual forward curve. 

%%%
% The first variable argument is the limit above which drain current should
% be ignored. This is due to a current-limited source being used causing
% the observed Id vs Vd behavior to no longer reflect the device behavior. 
id_limit = varargin{1}{:};

%%%
% Read data from the provided .csv file into a table and identify which
% columns correspond to gate voltage, drain voltage, and drain currents. 
raw_data = read_data(data_file);
headings = lower(raw_data(1,:));
vg_col = strcmpi(headings, 'v_g');
vd_col = strcmpi(headings, 'v_d');
id_matches = strfind(headings, 'i_d');
id_cols = ~cellfun('isempty', id_matches);

%%%
% Extract the gate voltage values, the drain voltage values, and the drain
% current values. The gate voltage vector corresponds to the to the fixed
% gate voltages at which forward curves were taken. The drain voltage
% vector corresonds to the values at which drain current was measured for
% each forward curve. The drain current matrix is organized such that each
% column corresponds to a unique forward curve. id and vg are organized 
% such that id(:,i) was measured at vg(i).
vg = raw_data{2, vg_col};
vd = raw_data{2, vd_col};
id = vertcat(raw_data{2, id_cols});

%%%
% Remove data above the current limit and report data in a cell array
% organized as {vg, vd, id}. 
id(id >= id_limit) = NaN;

data = cell(1,3);
data{1,1} = vg;
data{1,2} = vd;
data{1,3} = id;
end

