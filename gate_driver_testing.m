folder = fileparts(which(mfilename));
addpath(genpath(folder));
warning('off', 'all');

datafile = 'T4E3';
load(datafile)

figure
subplot(3, 1, 1)
plot(t_emp, vds_emp)
subplot(3, 1, 2)
plot(t_emp, id_emp)
subplot(3, 1, 3)
plot(t_emp, vgs_emp)