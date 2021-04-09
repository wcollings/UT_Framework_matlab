%%% Automated Model Generation using Static and Dynamic Data
% A model is autonomously optimized to predict static data (IV and CV) and
% dynamic data (DPT switching waveforms). User controls are intended to be
% input via configuration files. The 'config_master.txt' configuration file
% indicates which other configuration files should be used for sequential
% processing stages.


%%% Tidy Workspace
% Clear variales, close all figures, and clear command line. This is not
% required but is intented to reduce clutter. Warnings are turned off due
% to the use of UTF_16LE variables in LTspice2Matlab.m which is no longer a
% supported data type.
clear; close all; clc;
warning('off', 'all');
t_start = tic;


%%% Add nested folders to path
% Nested folder locations are collected and added to path in order to allow
% the current working folder to be more freely customized by the user. 
folder = fileparts(which(mfilename));
addpath(genpath(folder));

 
%%% Control optimization directly 
% Due to time constraints, new framework version will not utilize config
% files and will instead require user input and control directly. A new
% framework version will be optimized for realistic use cases. 
% par_names = {'Ld', 'Lg', 'Ls', 'kgd', 'kgs', 'kds'};
% par_vals = [0.15, 4.895, 2.6, .85, .85, .85];
% opt_vals = [31.2750    8.8950   11.1000    0.8500    1.1000    0.8500];
par_names = {'Ld', 'Lg', 'Ls', 'kgd', 'kgs', 'kds', 'Rd', 'Rg'};
par_vals = [0.15, 4.895, 2.6, .85, .85, .85, 5e-3, 1.5];
% par_vals = [0.1500    4.8950    2.6000    0.8500    0.8500    0.8500    0.0050    1.5000];
lb = [0.1 0.1 0.1 .7 .7 .7 1e-3 1e-3];
ub = [1000 1000 1000 2 2 2 1e3 1e3];

% cost = cost_handler(vals);
options = optimoptions('patternsearch','MaxTime',30*60,'MeshTolerance',1e-12);
optim_vals = patternsearch(@cost_handler, par_vals, [], [], [], [], lb, ub, [], options);
final_cost = cost_handler(optim_vals);


%%% Report framework duration
total_time = toc(t_start); 
disp(' ')
disp(['Total runtime: ', num2str(total_time)])