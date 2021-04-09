function [x, y] = objective_handler(test, event, signal)
%OBJECTIVE_HANDLER takes simulation file and returns simulation data in a
%form that is comparable to test data.

path = 'C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\SPICE_Files\';
file = ['Test', num2str(test)];
source = [path, file];
ltspicepath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
circuit = [source, '.asc'];
command = ['"',ltspicepath,'"', ' -run -b',' ','"',circuit,'"'];
% jsystem(command);
system(command);

rawdata = LTspice2Matlab(source);
x_raw = rawdata.time_vect; 

sigmap = {'vds', 'V(drain)'; 'vgs', 'V(gate)'; 'id', 'I(Lbus)'};
map_loc = strcmpi(sigmap(:,1), signal);
sig_name = sigmap{map_loc,2};
headings = rawdata.variable_name_list; 
sig_loc = strcmpi(headings, sig_name); 

data = rawdata.variable_mat; 
y_raw = data(sig_loc,:); 

% dxmin = min(diff(x_raw)); 
% xnew = x_raw(1):dxmin:x_raw(end); 
% ynew = pchip(x_raw, y_raw, xnew);

n = 10000;
xnew = linspace(x_raw(1), x_raw(end), n);
ynew = pchip(x_raw, y_raw, xnew);

cp = findchangepts(ynew, 'maxnumchanges', 4); 
eidx = fix(mean([cp(event-1:event);cp(event:event+1)]));
x1 = xnew(eidx(1)); 
x2 = xnew(eidx(2)); 
idx = x_raw>=x1&x_raw<=x2; 
x = x_raw(idx);
y = y_raw(idx);  

% cp = findchangepts(y_raw, 'maxnumchanges', 4); 
% eidx = fix(mean([cp(event-1:event);cp(event:event+1)]));
% x = x_raw(eidx(1):eidx(2));
% y = y_raw(eidx(1):eidx(2)); 


end

