function [opt_vals] = pattern(data,pars)
    %This quick and dirty version of PATTERN takes minimal imputs and contains
    %some hard-coded values in order to expedite paper completion. This will be
    %re-worked in a future version.
    %
    % function [outputArg1,outputArg2] = pattern(data,pars,obj,oarg,cost,carg,aarg)
    %PATTERN takes in a set of empirical data, a set of parameters, and uses a
    %pattern search algorithm to choose the parameter values which best predict
    %the empirical data in a SPICE simulation.
    %   data is a cell array with the first row containing headings and the
    %   second row containing column vectors of data. pars is a cell array with
    %   the first column containing parameter names and the second containing
    %   a vector of the initial parameter value, the percentage to sweep in the
    %   negative direction, and the percentage to sweep in the positive
    %   direction. obj is a handle for the objective function where oarg
    %   contains any arguments with which to call the objective function. cost
    %   is the cost function with carg containing any relevant cost arguments.
    %   aarg contains any relevant arguments to be passed to the algorithm.

    %%%%%%%%%%%%%
    % tran_data = isolate_transient(data,'vds',4);
    % sig_t = tran_data{3,1};
    % sig_d = tran_data{3,2};
    tran_data = isolate_transient(data,'id',3); %This was "4", but the data they gave us this time was incomplete
    sig_t = tran_data{3,1};
    sig_d = tran_data{3,2};
    global emp_pos emp_f emp_alpha par_names circuit model ltspicepath source
    [emp_pos, emp_f, emp_alpha] = template_pars(sig_t, sig_d);
    par_names = pars(:,1);
    %file = 'test1';
    %path = 'SPICE_Files\';
    %source = [path, file];
    [path,source,ext]=fileparts(circuit);
    source=strcat(path,"\",source,ext);
    %ltspicepath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
    %circuit = [source, '.asc'];
    %model = [path, 'Device_Models\Fraunhofer.lib'];
    %model = [path, 'Device_Models\EPC2022.lib'];
    %%%%%%%%%%%%%

    par_values = [pars{:,2}];
    par_values = reshape(par_values,[3,size(pars,1)]);

    static_values = par_values(1,:);
    min_values = static_values.*(1-par_values(2,:)/100);
    max_values = static_values.*(1+par_values(3,:)/100);
    %display=off or iter
    options = optimoptions('patternsearch','Display','iter', 'TolMesh', 1e-12,'TolX', 1e-20, 'TolCon', 1e-20,'PlotFcn',@psplotbestf);
    opt_vals = patternsearch(@quickSPICEcost,static_values,[],[],[],[],min_values,max_values,[],options);
    % options = optimset('PlotFcns',@optimplotfval);
    % opt_vals = fminsearch(@quickSPICEcost,static_values,options);

    % Write all parameter values to library file
    for k = 1:length(opt_vals)
        setpar(circuit,par_names{k},opt_vals(k));
    end

end

function cost = quickSPICEcost(parvals)
    
    global emp_pos emp_f emp_alpha par_names circuit model ltspicepath source
    
    % Write all parameter values to library file
    
    %this just prints the parameter values being tested
    s=sprintf("%s=%0.5e", par_names{1},parvals(1));
    for j=2:length(parvals)
        s=strcat(s,sprintf(", %s=%0.5e", par_names{j},parvals(j)));
    end
    %disp(s);
    
    for k = 1:length(parvals)
        setpar(circuit,par_names{k},parvals(k)); % should be "model", not "circuit"
    end
    
    % Run SPICE simulation
    command = ['"',ltspicepath,'"', ' -run -b ',circuit];
    jsystem(strjoin(command, ''));
    [path,fname,ext]=fileparts(source);
    rawf=strcat(path,'\',fname,'.raw');
    rawdata = LTspice2Matlab(rawf);
    headings = rawdata.variable_name_list;
    
    time_vect = rawdata.time_vect;
%     sim_idx = strcmpi(headings,'V(drain)');
%     vds_sim = rawdata.variable_mat(sim_idx,:);
%     dt_min = 0.4e-9; 
%     t_new = time_vect(1):dt_min:time_vect(end);
%     vds_samp = interp1(time_vect,vds_sim,t_new);
%     
%     sim_data = {'time','vds';t_new,vds_samp};
%     
%     tran_sim = isolate_transient(sim_data,'vds',4);
%     sim_t = tran_sim{3,1};
%     sim_d = tran_sim{3,2};
    sim_idx = strcmpi(headings,'Ix(m1:DRAININ)');
    id_sim = rawdata.variable_mat(sim_idx,:);
    dt_min = 0.4e-9; 
    t_new = time_vect(1):dt_min:time_vect(end);
    id_samp = interp1(time_vect,id_sim,t_new);
    
% % %     figure(2)
% % %     clf
% % %     plot(t_new, id_samp)
    
% % % % %     sim_data = {'time','vds';t_new,id_samp};
    sim_data = {'time','id';t_new,id_samp};
    
    tran_sim = isolate_transient(sim_data,'id',3);
%     sim_t = tran_sim{3,1};
%     sim_d = tran_sim{3,2};    
    sim_t = tran_sim{4,1}; %change these for different waveforms from simulation
    sim_d = tran_sim{4,2}; %change these for different waveforms from simulation
    
    figure(2)
    clf
    plot(sim_t, sim_d)
    
    [sim_pos, sim_f, sim_alpha] = template_pars(sim_t, sim_d);
    
    c_pos = abs((emp_pos-sim_pos)/emp_pos);
    c_f = abs((emp_f-sim_f)/emp_f);
    c_alpha = abs((emp_alpha-sim_alpha)/emp_alpha);
    
    cost = mean([c_pos, c_f, c_alpha]);
end

function setpar(netlist, parname, parval)
    fcell = {};
    fid = fopen(netlist,'r');
    nextline = fgetl(fid);
    n_line = 1;
    while nextline~=-1
        targidx = regexpi(nextline, parname);
        if ~isempty(targidx)
            curline=strsplit(nextline);
            curline{end}=num2str(parval);
            nextline=strjoin(curline);
        end
%         targidx = regexpi(nextline, ['.param ', parname]);
%         if ~isempty(targidx)
%             curline = strsplit(nextline);
%             par_loc = strcmpi(curline, '=');
%             targ_id = find(par_loc) + 1;
%             curline{targ_id} = num2str(parval);
%             nextline = strjoin(curline);
%        end
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
