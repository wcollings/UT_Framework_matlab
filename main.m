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


%%% Read master configuration file
% The master configuration file 'config_master.txt' will contain
% information regarding the stages of fitting. Once read in as a structure,
% each field corresponds to a fitting stage with a unique configuration
% file.
master_config_file = 'Configuration_Files/config_master.txt';
master_config = ini2struct(master_config_file);
stages = fieldnames(master_config);
n_stages = numel(stages);


%%% Perform each fitting stage serially
for s = 1:n_stages
    %%%
    % Extract information from master configuration file: 
    % The currient field corresponds to the current stage of fitting. The
    % field should contain the name of its configuration file and its data
    % source file.
    cur_stage = stages{s};
    cur_field = master_config.(cur_stage);
    cur_config_file = cur_field.config;
    cur_data_file = cur_field.data; 
    
    %%%
    % Extract information from current configuration file: 
    %   * Read the configuration file as a structure.
    %   * Extract modules for the objective and cost functions and the
    %     handling and output processes.
    %   * Extract arguments for the extracted modules.
    %   * Extract data from the specified source location
    %   * Extract the optimization algorithm to be used and any arguments.
    %   * Extract the fit parameters along with their search space data. 
    % In order to extract data from its source, the handling function is
    % used. This means that the 'data' variable can only be populated
    % following extraction of the handling function and its arguments. The
    % data is formatted by the handling function such that it is compatible
    % with the cost function being used. 
    config = ini2struct(cur_config_file);
    modules = config.modules;
    args = config.arguments;
    
    %%% Identify modules being used
    handling = eval(modules.handling);
    obj = eval(modules.objective);
    cost = eval(modules.cost);
    output = eval(modules.output);
    alg_params = config.algorithm_parameters;
    algorithm = eval(alg_params.algorithm);
    
    %%% Implement the specified handling process to acquire data in the
    %%% desired form
    handling_args = eval(args.handling);
    data = handling(cur_data_file, handling_args);
    
    %%% Identify the parameters to be used for optimization
    par_names = fieldnames(config.parameters_to_optimize);
    num_par = numel(par_names);
    pars = cell(num_par, 2);
    pars(:,1) = par_names;
    for p = 1:num_par
        cur_par_str = ['config.parameters_to_optimize.',par_names{p}];
        cur_par_rawdata = eval(cur_par_str);
        cur_par_data = strsplit(cur_par_rawdata);
        par_vals = cellfun(@str2num,cur_par_data);
        pars(p, 2) = {par_vals};
    end
    global circuit model ltspicepath source;
    %circuit=args.objective(2);
    %model=args.objective(3);
    %ltspicepath=args.objective(4);
    
    %%%
    % Implement the specified optimization algorithm.
    obj_args = eval(args.objective);
    circuit=string(obj_args(2));
    model=string(obj_args(3));
    ltspicepath=string(obj_args(4));
    cost_args = eval(args.cost);
    alg_args = eval(args.algorithm);
    %output_vars = algorithm(data, pars, obj, obj_args, cost, cost_args,... 
    %     alg_args);
    output_vars = algorithm(data, pars);
    % Output the relevant data as specified by the output process function.
    %output_args = eval(args.output);
    %output(output_vars, output_args); 
end


%%% Report framework duration
total_time = toc(t_start); 
disp(' ')
disp(['Total runtime: ', num2str(total_time)])
 par_vals=[ reshape([ pars{:,2} ], [3,size(pars,1)]) ];
 s=sprintf("%s=%0.5f ", par_names{1},par_vals(1,1));
    for j=2:length(par_vals(:,1))
      s=strcat(s,sprintf(", %s=%0.5e", par_names{j},par_vals(1,j)));
    end
  disp(s);

%%%%%%%%%%%%%%%%%%%%%
%Updates a parameter in the netlist file
%
%Takes the netlist itself,
%The parameter to change,
%And the value to set it to
%%%%%%%%%%%%%%%%%%%%%%
function setpar(netlist, parname, parval)
    
    fcell = {};
    fid = fopen(netlist,'r');
    nextline = fgetl(fid);
    n_line = 1;
    while nextline~=-1
        %targidx = regexpi(nextline, '.param VDS');
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
