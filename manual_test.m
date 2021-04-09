folder = fileparts(which(mfilename));
addpath(genpath(folder));
warning('off', 'all');

par_names = {'Ld', 'Lg', 'Ls', 'kgd', 'kgs', 'kds', 'Rd', 'Rg'};
% par_vals = [5, 10, 3, 1, 1, 1, 5e-3, 1];
% par_vals = [0.1500    4.8950    2.6000    0.8500    0.8500    0.8500    0.0050    1.5000];
% par_vals = [0.1500    4.8950    2.6000    1    1.5    1    0.0050    1.5000];
% par_vals = [0.1500    8    2.6000    .7    1    .7    0.0050    5];
par_vals = [2    6    2.6000    .7    1    .7    0.005    5];

% xrange = [2.51e-5 2.545e-5];
xrange = [3.035e-5 3.065e-5];


cost = cost_handler(par_vals);

plot_overlay

subplot(3, 1, 1)
xlim(xrange);
subplot(3, 1, 2)
xlim(xrange);
subplot(3, 1, 3)
xlim(xrange);