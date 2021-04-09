folder = fileparts(which(mfilename));
addpath(genpath(folder));


% test = 4;
% event = 2;
% t_shift = -.382e-5;

test = 4; 
event = 3; 
% t_shift = -.382e-5; 
t_shift = -.38e-5; 

load(['T',num2str(test),'E',num2str(event)])

path = 'C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\SPICE_Files\';
file = ['Test', num2str(test)];
source = [path, file];
ltspicepath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
circuit = [source, '.asc'];
command = ['"',ltspicepath,'"', ' -run -b',' ',circuit];
% jsystem(command);
system(command);

rawdata = LTspice2Matlab(source);
t_sim = rawdata.time_vect; 

% sigmap = {'vds', 'V(drain)'; 'vgs', 'V(gate)'; 'id', 'I(Lbus)'};
% map_loc = strcmpi(sigmap(:,1), signal);
% sig_name = sigmap{map_loc,2};
headings = rawdata.variable_name_list; 
vds_loc = strcmpi(headings, 'V(drain)'); 
vgs_loc = strcmpi(headings, 'V(gate)'); 
id_loc = strcmpi(headings, 'I(Lbus)'); 

data = rawdata.variable_mat; 
vds_sim = data(vds_loc,:); 
vgs_sim = data(vgs_loc,:);
id_sim = data(id_loc,:);
% id_sim = data(id_loc,:)-3;

lw = 1.2; 
figure
subplot(3, 1, 1)
hold on
plot(t_emp, vds_emp, '-k', 'linewidth', lw)
plot(t_sim-t_shift, vds_sim, ':r', 'linewidth', lw)
hold off

subplot(3, 1, 2)
hold on
plot(t_emp, id_emp, '-k', 'linewidth', lw)
plot(t_sim-t_shift, id_sim, ':r', 'linewidth', lw)
hold off

subplot(3, 1, 3)
hold on
plot(t_emp, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim-t_shift, vgs_sim, ':r', 'linewidth', lw)
hold off

