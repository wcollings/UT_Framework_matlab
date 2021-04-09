folder = fileparts(which(mfilename));
addpath(genpath(folder));
warning('off','all')


align_by_event = 2;


handling_args = {11};
data = handling_transient('Test3.csv', handling_args);
t_emp = data{2,1};
vds_emp = data{2,2};
vgs_emp = data{2,3};
id_emp = data{2,4};


event2 = [20.1 20.4];
event3 = [25.3 25.6];



%%%%%%%%%% Optimized values
par_names =     {'ld','r_dp','c_dp','lg','r_gp','c_gp','ls','r_sp','c_sp','rd','rg','rs','kds','kgs','kdg','r_ds','r_gs','r_dg'};
par_vals_stat = [0.1,1e6,0,4.895,1e6,0,2.6,1e6,0,1e-4,1e-4,1e-4,1.0,1.0,1.0,1e-4,1e-4,1e-4];
%%%%%%%%% optimization_costvalues_2020_02_05_1652
par_vals_dyn = [0.10391 18 932.875 5.895 9.749 389 3.6 4 401.0005 0.0001 0.0001 0.0001 1 3.002 1 1e-4 1e-4 1e-4];




% Define file paths
file = 'test3';
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



%%%%%%%%%% Event 2
manual_shift_stat2 = 0;
manual_shift_stat3 = -0.009;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2, 3, 1)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2+manual_shift_stat2, vds_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim(event2)
% ylim([-5 800])
set(gca, 'fontweight', 'bold')
% 
subplot(2, 3, 2)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2+manual_shift_stat2, id_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim(event2)
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 3)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static2+manual_shift_stat2, vgs_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim(event2)
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
xlim(event2)
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
xlim(event2)
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
xlim(event2)
% ylim([-5 800])
set(gca, 'fontweight', 'bold')


%%%%%%%%%% Event 3
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2, 3, 1)
hold on
plot(t_emp*1e6, vds_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3+manual_shift_stat3, vds_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_D_S (V)')
xlim(event3)
ylim([-5 800])
set(gca, 'fontweight', 'bold')
% 
subplot(2, 3, 2)
hold on
plot(t_emp*1e6, id_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3+manual_shift_stat3, id_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('I_D (A)')
xlim(event3)
% ylim([-5 800])
set(gca, 'fontweight', 'bold')

subplot(2, 3, 3)
hold on
plot(t_emp*1e6, vgs_emp, '-k', 'linewidth', lw)
plot(t_sim_stat*1e6+t_shift_static3+manual_shift_stat3, vgs_sim_stat, ':r', 'linewidth', lw)
hold off
grid on
box on
xlabel('Time (\mus)')
ylabel('V_G_S (V)')
xlim(event3)
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
xlim(event3)
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
xlim(event3)
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
xlim(event3)
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