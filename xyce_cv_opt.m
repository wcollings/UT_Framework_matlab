%%% Read data from a .prn function of the form output by Xyce
% xyce_path = 'C:\Program Files\Xyce 6.11.1 OPENSOURCE\bin\Xyce.exe';

folder = fileparts(which(mfilename));
addpath(genpath(folder));

% pars_gd = {'cogd', 'xgd1', 'xgd2', 'agd', 'chigd'};
vals_gd = [1.4540e-09 1.5041e+04 0.0025 1 .5];
lb = [1e-12 1 0 0 0.1];
ub = [1e-6 1e6 1 10 0.9];

cost1 = cost_gd(vals_gd);

opt_vals = patternsearch(@cost_gd, vals_gd,[],[],[],[],lb,ub);
% opt_vals = fminsearch(@cost_gd, vals_gd);

cost2 = cost_gd(opt_vals);

function cost = cost_gd(vals_gd)
    pars_gd = {'cogd', 'xgd1', 'xgd2', 'agd', 'chigd'};
    netlist = 'C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\SPICE_Files\SICFET.lib';
    for k = 1:numel(pars_gd)
        update_netlist(netlist, pars_gd{k}, vals_gd(k));
    end
    c = get_xyce_cv();
    load('C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\Fitting_Data\Static_Data\cv_pico_logsamp.mat')
    vds_start = 0;
    vds_end = log10(vds_emp(end));
    n_pts = 30;
    cgd_sim = c(3,:); 
    cgd_emp = c_emp(3,:); 
%     vds_sweepvec = logspace(vds_start, vds_end, n_pts);
%     vds_sweepvec = logspace(vds_end, vds_start, n_pts);
    vds_sweepvec = linspace(10.^vds_start, 10.^vds_end, n_pts);
%     c_interp = pchip(vds_sweepvec, c, vds_emp);
    c_interp = pchip(vds_emp, cgd_emp, vds_sweepvec);
    cost = norm(cgd_sim-c_interp);
end

function update_netlist(netlist, parameter, val)
    fcell = {};
    fid = fopen(netlist,'r');
    nextline = fgetl(fid);
    n_line = 1; 
    targstr = ['.param ', parameter];
    while nextline ~= -1
        targidx = regexpi(nextline, targstr);
        if ~isempty(targidx)
            curline = strsplit(nextline);
            vds_loc = strcmpi(curline, parameter);
            targ_id = find(vds_loc) + 2; 
            curline{targ_id} = num2str(val);
            nextline = strjoin(curline);
        end
        fcell{n_line} = nextline; 
        n_line = n_line + 1;
        nextline = fgetl(fid);
        if isempty(nextline)
            nextline = ' ';
        end
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
