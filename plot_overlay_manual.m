folder = fileparts(which(mfilename));
addpath(genpath(folder));
warning('off','all')

% raw_data = readtable('test1.csv');

align_by_event = 2;

% Collect empirical data
% handling_args = {11,{'time','vds'}};
handling_args = {11};
data = handling_transient('Test1.csv', handling_args);
t_emp = data{2,1};
vds_emp = data{2,2};
vgs_emp = data{2,3};
id_emp = data{2,4};

% Define model tuning parameter values
% par_names = {'ld','lg','ls','kds','kgs','kdg'};
% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1];
% % par_vals_dyn = [0.099985, 4.6294, 1, 0.70313, 1];
% par_vals_dyn = [0.1, 4.895, 2.6, 1, 3.6, 1];
% % .PARAM Ld = 0.099985
% % .PARAM Lg = 4.6294
% % .PARAM Ls = 3.1076
% % .PARAM kgs = 1
% % .PARAM kds = 0.70313
% % .PARAM kdg = 1
% % % .PARAM Ld = 0.1
% % % .PARAM Lg = 5.4712
% % % .PARAM Ls = 3.0687
% % % .PARAM kgs = 1
% % % .PARAM kds = 0.74609
% % % .PARAM kdg = 0.75787
% par_names = {'ld','lg','ls','kds','kgs','kdg','rd','rg','rs','r_dp','r_ds','r_gs'};
% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1, 1, 1, 1e-3, 10];
% par_vals_dyn = [0.1, 4.7075, 2.6, 1, 1, 1, 2.9375, 7, .001, 10.0039];
% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10];
% par_vals_dyn = [0.1, 4.7075, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10.0039];
% par_vals_dyn = [0.096094, 3.977, 3.1078, 0.99994, 0.71094, 0.95313, 1e-3, 1e-3, 1e-3, 9.0625];
% par_vals_dyn = [12.125, 14, 5.25, 1, .875, 1, 1e-3, 1e-3, 1e-3, 11.5];
% par_vals_dyn = [0.1, 4.895, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10];

% par_vals_dyn = [0.1, 4.8639, 2.6, 1, 1, 1, 1.25, 2.25, 0.001, 16];
% par_vals_dyn = [0.098047 3.082 3.85 1 0.99609 1 1.2188 10.5 0.001 5];
% par_vals_dyn = [5 3.082 3.85 1 0.99609 1 2 10.5 0.001 1e6];
% par_vals_dyn = [0.1 3.645 3.85 0.875 1 1 1.1328 4.5 0.001 16.6406];


% %%%%%%%%%%%%
% par_vals_dyn = [0.1 3.645 3.85 0.875 5 1 1.1328 22 0.001 16.6406];
% %%%%%%%%%%%%

% %%%%%%%%%%%%
% par_vals_dyn = [0.099939 3.895 4.1 0.9375 5.375 1 1.0395 22 0.001 15.2656];
% %%%%%%%%%%%%

% %%%% optimization_costvalues_2020_01_30_1049
% par_vals_dyn = [0.10098 4.895 3.6 1 6 0.5 1 5 0.001 46];
% %%%%

% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10, 1];
% %%%%%%%%%%% optimization_costvalues_2020_01_30_1643
% par_vals_dyn = [0.1 4.6606 3.7245 0.9375 6 1 0.5 5.501 0.001 5.0078 0.49805, 1e-4];
% %%%%%%%%%%%

% %%%%%%%%%%%% optimization_costvalues_2020_01_31_1620
% par_vals_dyn = [0.10391 4.645 3.849 1.002 2.5001 1 0.5 5.501 0.001 98.3086 0.375 10.749];
% %%%%%%%%%%%%

% %%%%%%%%%%%%
% par_vals_dyn = [0.1 4.6606 3.7245 0.9375 1 1 0.5 5.501 0.001 5.0078 0.49805];
% %%%%%%%%%%%%
% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10, 1, 1];
% par_vals_dyn = [0.1 3.645 3.8969 0.99219 1.0938 1.0625 0.5 5.501 0.001 37.75 0.90625 10.875];


%%%% Manual Tuning
% par_names =     {'ld',  'lg',   'ls',   'kds',  'kgs',  'kdg',  'rd',   'rg',   'rs',   'r_dp', 'r_ds', 'r_gs', 'r_dg', 'c_dp'};
% par_vals_stat = [0.1,   4.895,  2.6,    1,      1,      1,      1e-3,   1e-3,   1e-3,   1e6,    1e-4,   1e-4,   1e-4,   0];
% par_vals_dyn =  [0.1,   4.895,  2.6,    1,      1,      1,      1e-3,   1e-3,   1e-3,   1,      1e-6,   1e-6,   1e-4,   1e-9];

% par_names =     {'ld',  'r_dp', 'c_dp', 'lg',   'r_gp', 'c_gp', 'ls',   'r_sp', 'c_sp', 'rd',   'rg',   'rs',   'kds',  'kgs',  'kdg',  'r_ds', 'r_gs', 'r_dg'};
% par_vals_stat = [0.1,   1e6,    0,      4.895,  1e6,    0,      2.6,    1e6,    0,      1e-4,   1e-4,   1e-4,   1.0,    1.0,    1.0,    1e-4,   1e-4,   1e-4];
% par_vals_dyn =  [0.1,   10,    1,   4.895,  10,     1,   2.6,    20,     1,   1e-4,   1e-4,   1e-4,   1.0,    1.0,    1.0,    1e-4,   1e-4,   1e-4];
% % 
% %%%% optimization_costvalues_2020_02_03_1227
% par_vals_dyn = [1.1 17.998 1 6.7705 2.0508 1 2.5922 2 1 0.0001 0.0001 0.0001 .9091 0.9830 .9524 1.9991 0.03135 0.0001];

%%%%%%%%%%%% optimization_costvalues_2020_01_31_1620
% par_names = {'ld','lg','ls','kds','kgs','kdg','rd','rg','rs','r_dp','r_ds','r_gs'};
% par_vals_dyn = [0.10391 4.645 3.849 1.002 2.5001 1 0.5 5.501 0.001 98.3086 0.375 10.749];
% % % par_names =     {'ld','r_dp','c_dp','lg','r_gp','c_gp','ls','r_sp','c_sp','rd','rg','rs','kds','kgs','kdg','r_ds','r_gs','r_dg'};
% % % par_names = {'ld','lg','ls','kds','kgs','kdg','rd','rg','rs','r_dp','r_ds','r_gs'};
% % % par_vals_dyn = [0.10391 () () 4.645 () () 3.849 () () () () () 1.002 2.5001 1 0.5 5.501 0.001 98.3086 0.375 10.749];
%%%%%%%%%%%%

% % %%%%
% %%%%% optimization_costvalues_2020_02_04_0524
% par_vals_dyn = [0.11563 20 945 6.895 18 1000.9961 3.6 4 977 0.0001 0.0001 0.0001 1 1 0.99998 0.0001 8.0001 0.0001];

% .PARAM Ld = 0.096094
% .PARAM Lg = 3.977
% .PARAM Ls = 3.1078
% .PARAM kgs = 0.99994
% .PARAM kds = 0.71094
% .PARAM kdg = 0.95313
% .PARAM rdp = 9.0625

% %%%%% Best Event 3
% %%%%%%%%%% optimization_costvalues_2020_02_04_1721
% par_vals_dyn = [0.225 9 954 2.895 20 144.25 3.85 8785.25 1000.5 0.0001 0.0001 0.0001 0.9375 3.0625 1 0.7501 4.5001 0.0001];
% %%%%%%%%%%


% .PARAM Ld = 0.225
% .PARAM r_dp = 9
% .PARAM c_dp = 954
%  
% .PARAM Lg = 2.895
% .PARAM r_gp = 20
% .PARAM c_gp = 144.25
%  
% .PARAM Ls = 3.85
% .PARAM r_sp = 8785.25
% .PARAM c_sp = 1000.5
%  
% .PARAM rd = 0.0001
% .PARAM rg = 0.0001
% .PARAM rs = 0.0001
%  
% .PARAM kds = 0.9375
% .PARAM kgs = 3.0625
% .PARAM kdg = 1
% .PARAM r_ds = 0.7501
% .PARAM r_gs = 4.5001
% .PARAM r_dg = 0.0001


% %%%%%%%% optimization_costvalues_2020_02_05_1324 (optimized for id 3)
% par_vals_dyn = [0.10391 20 1.4922 4.895 10 904.75 3.5922 26774.625 961.9844 0.0001 0.0001 0.0001 1 1 1 7.5626 0.0001 0.0001];
% %%%%%%%%
% .PARAM Ld = 0.10391
% .PARAM r_dp = 20
% .PARAM c_dp = 1.4922
%  
% .PARAM Lg = 4.895
% .PARAM r_gp = 10
% .PARAM c_gp = 904.75
%  
% .PARAM Ls = 3.5922
% .PARAM r_sp = 26774.625
% .PARAM c_sp = 961.9844
%  
% .PARAM rd = 0.0001
% .PARAM rg = 0.0001
% .PARAM rs = 0.0001
%  
% .PARAM kds = 1
% .PARAM kgs = 1
% .PARAM kdg = 1
% .PARAM r_ds = 7.5626
% .PARAM r_gs = 0.0001
% .PARAM r_dg = 0.0001


% %%%%%%%%% optimization_costvalues_2020_02_05_1652
% par_vals_dyn = [0.10391 18 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3.002 1 0.0001 0.0001 0.0001];
% %%%%%%%%%%


% .PARAM Ld = 0.10391
% .PARAM r_dp = 18
% .PARAM c_dp = 932.875
%  
% .PARAM Lg = 5.895
% .PARAM r_gp = 9.749
% .PARAM c_gp = 389
%  
% .PARAM Ls = 3.6
% .PARAM r_sp = 4
% .PARAM c_sp = 401.0005
%  
% .PARAM rd = 0.0001
% .PARAM rg = 0.0001
% .PARAM rs = 0.0001
%  
% .PARAM kds = 1
% .PARAM kgs = 3.002
% .PARAM kdg = 1
% .PARAM r_ds = 0.0001
% .PARAM r_gs = 0.0001
% .PARAM r_dg = 0.0001








%%%%%%%%%% Optimized for turn-off Vds
% par_names = {'ld','lg','ls','kds','kgs','kdg','rd','rg','rs','r_dp','r_ds','r_gs'};
% par_vals_stat = [0.1, 4.895, 2.6, 1, 1, 1, 1e-3, 1e-3, 1e-3, 10, 1];
%%%%%%%%%%% optimization_costvalues_2020_01_30_1643
% par_vals_dyn = [0.1 4.6606 3.7245 0.9375 6 1 0.5 5.501 0.001 5.0078 0.49805, 1e-4];
%%%%%%%%%%%

% %%%%%%%%%% Optimized for turn-on Id
% par_names =     {'ld','r_dp','c_dp','lg','r_gp','c_gp','ls','r_sp','c_sp','rd','rg','rs','kds','kgs','kdg','r_ds','r_gs','r_dg'};
% par_vals_stat = [0.1,1e6,0,4.895,1e6,0,2.6,1e6,0,1e-4,1e-4,1e-4,1.0,1.0,1.0,1e-4,1e-4,1e-4];
% %%%%%%%%% optimization_costvalues_2020_02_05_1652
% par_vals_dyn = [0.10391 18 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3.002 1 0.0001 0.0001 0.0001];
% %%%%%%%%%%





%%%%%%%%%% Optimized values
par_names =     {'ld','r_dp','c_dp','lg','r_gp','c_gp','ls','r_sp','c_sp','rd','rg','rs','kds','kgs','kdg','r_ds','r_gs','r_dg'};
par_vals_stat = [0.1,1e6,0,4.895,1e6,0,2.6,1e6,0,1e-4,1e-4,1e-4,1.0,1.0,1.0,1e-4,1e-4,1e-4];
%%%%%%%%% optimization_costvalues_2020_02_05_1652
par_vals_dyn = [0.10391 18 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3.002 1 1e-4 1e-4 1e-4];
%%%%%%%%%%% optimization_costvalues_2020_01_30_1643
% par_vals_dyn = [0.1 5.0078 932.875 4.6606 9.749 389 3.7245 4 401.0005 0.5 5.501 0.001  0.9375 6 1 0.49805, 1e-4 0.0001];
%%%%%%%%%%

% par_vals_on = [0.10391 18 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3.002 1 0.0001 0.0001 0.0001];
% par_vals_off = [0.1 5.0078 932.875 4.6606 9.749 389 3.7245 4 401.0005 0.5 5.501 0.001  0.9375 6 1 0.49805, 1e-4 0.0001];
% par_vals = [par_vals_on;par_vals_off];
% % par_vals_dyn = mean(par_vals);
% % par_vals_dyn = [0.10391 5.0078 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3 1 0.0001 0.0001 0.0001];
% w_on = 2; 
% w_off = 1; 
% par_vals_dyn = (w_on*par_vals_on+w_off*par_vals_off)/(w_on+w_off);





% Define file paths
file = 'test1';
path = 'C:\Users\andre\Documents\Research\Automated_SPICE_Modeling\UT_Framework_3.0.1\SPICE_Files\';
source = [path, file];
ltspicepath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
circuit = [source, '.asc'];
model = [path, 'Device_Models\Fraunhofer.lib'];


% Write static values to library file
for k = 1:length(par_vals_stat)
    setpar(model,par_names{k},par_vals_stat(k));
end
% Run SPICE simulation
command = ['"',ltspicepath,'"', ' -run -b',' ',circuit];
jsystem(command);
rawdata = LTspice2Matlab(source);
headings = rawdata.variable_name_list;
t_sim_stat = rawdata.time_vect;
sim_vds_idx = strcmpi(headings,'V(drain)');
vds_sim_stat = rawdata.variable_mat(sim_vds_idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_vgs_idx = strcmpi(headings,'V(gate)');
vgs_sim_stat = rawdata.variable_mat(sim_vgs_idx,:);
sim_id_idx = strcmpi(headings,'I(Rbus)');
id_sim_stat = rawdata.variable_mat(sim_id_idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write dynamic values to library file
for k = 1:length(par_vals_dyn)
    setpar(model,par_names{k},par_vals_dyn(k));
end
% Run SPICE simulation
command = ['"',ltspicepath,'"', ' -run -b',' ',circuit];
jsystem(command);
rawdata = LTspice2Matlab(source);
headings = rawdata.variable_name_list;
t_sim_dyn = rawdata.time_vect;
sim_vds_idx = strcmpi(headings,'V(drain)');
vds_sim_dyn = rawdata.variable_mat(sim_vds_idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_vgs_idx = strcmpi(headings,'V(gate)');
vgs_sim_dyn = rawdata.variable_mat(sim_vgs_idx,:);
sim_id_idx = strcmpi(headings,'I(Rbus)');
id_sim_dyn = rawdata.variable_mat(sim_id_idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t_shift_static = 4.21;

changes_emp = findchangepts(vds_emp, 'maxnumchanges', 4);
changes_sta = findchangepts(vds_sim_stat, 'maxnumchanges', 4);
changes_dyn = findchangepts(vds_sim_dyn, 'maxnumchanges', 4);

t2_emp = t_emp(changes_emp(align_by_event))*1e6;
t2_sta = t_sim_stat(changes_sta(align_by_event))*1e6;
t2_dyn = t_sim_dyn(changes_dyn(align_by_event))*1e6;

t3_emp = t_emp(changes_emp(3))*1e6;
t3_sta = t_sim_stat(changes_sta(3))*1e6;
t3_dyn = t_sim_dyn(changes_dyn(3))*1e6;

t_shift_static2 = t2_emp - t2_sta; 
t_shift_dynamic2 = t2_emp - t2_dyn;
t_shift_static3 = t3_emp - t3_sta; 
t_shift_dynamic3 = t3_emp - t3_dyn;
lw = 1.2;

% figure
% hold on
% plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
% plot(t_sim_stat*1e6+t_shift_static, vds_sim_stat, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('V_D_S (V)')
% xlim([15.9 16.15])
% ylim([-5 800])
% set(gca, 'fontweight', 'bold')
% % 
% figure
% hold on
% plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
% plot(t_sim_stat*1e6+t_shift_static, id_sim_stat, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('I_D (A)')
% xlim([15.9 16.15])
% % ylim([-5 800])
% set(gca, 'fontweight', 'bold')
% 
% figure
% hold on
% plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
% plot(t_sim_stat*1e6+t_shift_static, vgs_sim_stat, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('V_G_S (V)')
% xlim([15.9 16.15])
% % ylim([-5 800])
% set(gca, 'fontweight', 'bold')
% 
% % t_shift_dynamic = 4.21;
% 
% figure
% hold on
% plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
% plot(t_sim_dyn*1e6+t_shift_dynamic, vds_sim_dyn, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('V_D_S (V)')
% xlim([15.9 16.15])
% ylim([-5 800])
% set(gca, 'fontweight', 'bold')
% 
% figure
% hold on
% plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
% plot(t_sim_dyn*1e6+t_shift_dynamic, id_sim_dyn, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('I_D (A)')
% xlim([15.9 16.15])
% % ylim([-5 800])
% set(gca, 'fontweight', 'bold')
% 
% figure
% hold on
% plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
% plot(t_sim_dyn*1e6+t_shift_dynamic, vgs_sim_dyn, ':r', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('Time (\mus)')
% ylabel('V_G_S (V)')
% xlim([15.9 16.15])
% % ylim([-5 800])
% set(gca, 'fontweight', 'bold')

%%%%%%%%%% Event 2
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2, 3, 1)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2, vds_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim([15.9 16.15])
ylim([-5 800])
set(gca, 'fontweight', 'bold')
% 
subplot(2, 3, 2)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2, id_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim([15.9 16.15])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 3)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2, vgs_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim([15.9 16.15])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

% t_shift_dynamic = 4.21;

subplot(2, 3, 4)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic2, vds_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim([15.9 16.15])
ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 5)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic2, id_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim([15.9 16.15])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 6)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic2, vgs_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim([15.9 16.15])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')


%%%%%%%%%% Event 3
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2, 3, 1)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3, vds_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim([21.1 21.4])
ylim([-5 800])
set(gca, 'fontweight', 'bold')
% 
subplot(2, 3, 2)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3, id_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim([21.1 21.4])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 3)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3, vgs_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim([21.1 21.4])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

% t_shift_dynamic = 4.21;

subplot(2, 3, 4)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic3, vds_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim([21.1 21.4])
ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 5)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic3, id_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim([21.1 21.4])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 6)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_dyn*1e6+t_shift_dynamic3, vgs_sim_dyn, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim([21.1 21.4])
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

function setpar(netlist, parname, parval)
    fcell = {};
    fid = fopen(netlist,'r');
    nextline = fgetl(fid);
    n_line = 1;
    while nextline~=-1
%         targidx = regexpi(nextline, '.param VDS');
        targidx = regexpi(nextline, ['.param ', parname]);
        if ~isempty(targidx)
            curline = strsplit(nextline);
            par_loc = strcmpi(curline, '=');
            targ_id = find(par_loc) + 1;
            curline{targ_id} = num2str(parval);
            nextline = strjoin(curline);
        end
        fcell{n_line} = nextline;
        n_line = n_line + 1;
        nextline = fgetl(fid);
    end
    fclose(fid);
    fid = fopen(netlist, 'w');
    for i = 1:numel(fcell)
        nextline = fcell{i};
        slashes = regexp(nextline, '\');
        if isempty(slashes)
            fprintf(fid, nextline);
        else
            n_slash = length(slashes);
            apstr = '';
            for k = 1:n_slash
                sidx = slashes(n_slash + 1 - k);
                str1 = nextline(1:sidx-1);
                str2 = nextline(sidx+1:end);
                nextline = [str1, '%s', str2];
                apstr = [apstr, ', ', '''', '\', ''''];
            end
            command = ['fprintf(fid, nextline',apstr,');'];
            eval(command);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end