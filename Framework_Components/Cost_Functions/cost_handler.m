function cost = cost_handler(vals)
%%% COST_HANDLER takes a given search node and validation data and returns
% a cost to the calling function. 
% 
test = 4; 
event = 2; 
signal = 'vds'; 

% test = 4; 
% event = 3; 
% signal = 'id'; 

%%% Pars here must match those in main_noconfig.m
% pars = {'Ld', 'Lg', 'Ls', 'kgd', 'kgs', 'kds'};
pars = {'Ld', 'Lg', 'Ls', 'kgd', 'kgs', 'kds', 'Rd', 'Rg'};

datafile = ['T',num2str(test),'E',num2str(event)];
load(datafile)

%%% Update the model
libfile = 'SPICE_Files/MSA12080A.lib';
% libfile = 'SPICE_Files/SICFET.lib';
for k = 1:numel(pars)
    cur_par = pars{k}; 
    cur_val = vals(k);
    update_lib(libfile, cur_par, cur_val); 
end

[x_sim, y_sim] = objective_handler(test, event, signal);

cp_emp = findchangepts(eval([signal, '_emp']), 'maxnumchanges', 1);
cp_sim = findchangepts(y_sim, 'maxnumchanges', 1);

t_cp = t_emp(cp_emp); 
x_cp = x_sim(cp_sim); 

t_emp_ref = t_emp - t_cp; 
x_sim_ref = x_sim - x_cp;

y_samp = pchip(x_sim_ref, y_sim, t_emp_ref);

% figure
% hold on
% plot(t_emp, vds_emp)
% plot(t_emp, y_samp)
% hold off

% cost = norm(vds_emp - y_samp);
cost = template_matching(t_emp, eval([signal, '_emp']), y_samp);
cost = cost*100; 
end