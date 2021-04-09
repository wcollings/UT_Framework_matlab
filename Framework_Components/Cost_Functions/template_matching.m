function [ cost ] = template_matching( time, emp, sim )
%TEMPLATE_MATCHING takes in two sets of data and reports back an associated
%cost by fitting each data with parameters describing 'second-order-like'
%behaviors such as the percent overshoot, the attenuation constant, and the
%ringing frequency. 
cost = 0;

pos_w = 1;
alpha_w = 1;
f_w = 1;

total_weight = sum([pos_w, alpha_w, f_w]);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract data %%%%%%%%%%%%%%%%%%%%%%%%%%
% meas_t = data{i,1};    % Time vector in first column
% meas_d = data{i,j};    % Data vector in jth column
% sim_t = simulation{i,1};      % Time vector in first column
% sim_d = simulation{i,j};      % Data vector in jth column
[m_pos,m_f,m_alpha] = template_pars(time,emp); % Template parameters for measured data
[s_pos,s_f,s_alpha] = template_pars(time,sim);   % Template parameters for simulated data
%         [m_pos,m_f,m_alpha,m_s] = template_pars(meas_t,meas_d); % Template parameters for measured data
%         [s_pos,s_f,s_alpha,s_s] = template_pars(sim_t,sim_d);   % Template parameters for simulated data
pos_err = abs(s_pos - m_pos)/abs(m_pos)*pos_w;      % Normalize and weight pos error

f_err = abs(s_f - m_f)/abs(m_f)*f_w;                % Normalize and weight frequency cost
%         s_spec = s_f(1,:);
%         s_freq = s_f(2,:);
%         m_spec = m_f(1,:);
%         m_freq = m_f(2,:);
%         [m_s,i_s] = max(s_spec);
%         [m_m,i_m] = max(m_spec);
%         s_df = s_freq(i_s);
%         m_df = m_freq(i_m);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%% Dominant Frequency %%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         f_err = abs(s_df - m_df)/abs(m_df)*f_w;                % Normalize and weight frequency cost
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Dominant Frequency and Magnitude %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         lambda = 1;
%         f_err = f_w*(lambda*(1-s_df/m_df).^2+(1-m_s/m_m).^2);         % Normalize weight and frequency cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha_err = abs(s_alpha - m_alpha)/abs(m_alpha)*alpha_w;    % Normalize and weight alpha error
%         s_err = abs(s_s - m_s)/abs(m_s)*s_w;
%         cur_err = pos_err + f_err + alpha_err + s_err;  % Sum errors
cur_err = pos_err + f_err + alpha_err;  % Sum errors

cur_err_norm = cur_err/total_weight;           % Normalize to total weight
cost = cost + cur_err_norm;           % Add to total error sum


end

