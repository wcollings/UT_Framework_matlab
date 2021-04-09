addpath("Interpreters")
ltspicepath="C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe";
circuit="SPICE_Files\epc2022_dpt.net";
command = ['"',ltspicepath,'"', ' -run -b ',circuit];
jsystem(strjoin(command, ''));

d=readmatrix("Fitting_Data\Transient_Data\Test1_trimmed.csv"); 
raw = LTspice2Matlab("SPICE_Files\epc2022_dpt.raw");
vgs_ind = strcmpi(raw.variable_name_list,'V(g)');
vds_ind = strcmpi(raw.variable_name_list,'V(d)');
id_ind  = strcmpi(raw.variable_name_list,'Ix(m1:DRAININ)');
s=[raw.time_vect ; raw.variable_mat(vgs_ind,:); raw.variable_mat(vds_ind,:); raw.variable_mat(id_ind,:)]';
plot(d(:,1),d(:,5))
%hold
%plot(s(:,1),s(:,4))
legend(["exp","sim"])
