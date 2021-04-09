%%% Read data from a .prn function of the form output by Xyce
% datasource = 'C:\Users\andre\Documents\Rick Work\SaberRD_demo\Datasheets\SiC_MOSFET\rohmsic.mat';
xyce_path = 'C:\Program Files\Xyce 6.11.1 OPENSOURCE\bin\Xyce.exe';
% circuit = 'C:\Users\andre\Documents\Typhoon\Xyce\sic_fwd.cir';
% circuit = 'C:\Users\andre\Documents\Typhoon\Xyce\angelov.cir';

folder = fileparts(which(mfilename));
addpath(genpath(folder));

%%%%%%%%%%%%%%%%%%%%%%%%%% Custom Model %%%%%%%%%%%%%%%%%%%%%%%%%%
% datasource = 'C:\Users\andre\Documents\Research\SiC_Modeling\CREE Model\cvdata.mat';
% load(datasource)

datasource = 'NewData_CVData_pf.csv';
rawcv = read_data(datasource); 
headings = rawcv(1,:); 
vds_loc = strcmpi(headings, 'vds'); 
ciss_loc = strcmpi(headings, 'ciss'); 
coss_loc = strcmpi(headings, 'coss');
crss_loc = strcmpi(headings, 'crss'); 
vds = rawcv{2,vds_loc};
ciss = rawcv{2,ciss_loc};
coss = rawcv{2,coss_loc};
crss = rawcv{2,crss_loc};
c = [ciss; coss; crss];
sigs = {'C_I_S_S', 'C_O_S_S', 'C_R_S_S'};

% circuit = 'C:\Users\andre\Documents\Typhoon\Xyce\cv_testbed.cir';
circuit = 'C:\Users\andre\Documents\Research\Automated SPICE Modeling\cv_testbed.cir';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% RSiC Model %%%%%%%%%%%%%%%%%%%%%%%%%%%
% datasource = 'C:\Users\andre\Documents\Rick Work\SaberRD_demo\Datasheets\SiC_MOSFET\cv.mat';
% circuit = 'C:\Users\andre\Documents\Typhoon\Xyce\cv_testbed_rsic.cir';
% % Load datasheet data
% load(datasource)
% c = [ciss;coss;crss];
% c = c*1e12;
% sigs = {'C_I_S_S','C_O_S_S','C_R_S_S'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vds_emp = vds;
c_emp = c;
sigs_emp = sigs;
ciss_idx = strcmpi(sigs, 'C_I_S_S');
coss_idx = strcmpi(sigs, 'C_O_S_S');
crss_idx = strcmpi(sigs, 'C_R_S_S');
ciss_emp = c_emp(ciss_idx,:);
coss_emp = c_emp(coss_idx,:);
crss_emp = c_emp(crss_idx,:);
cgd_emp = crss_emp;
cds_emp = coss_emp - crss_emp;
cgs_emp = ciss_emp - crss_emp;

%%%%%%%%%%%%%%%%%%%%
filename = [circuit, '.FD.prn'];
delete(filename)
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vds_start = 0;
vds_end = log10(vds(end));
n_pts = 30;
vds_sweepvec = logspace(vds_start, vds_end, n_pts);
lw = 1.2;
fsz = 10;
figure(1)
clf
hold on
plot(vds_emp, ciss_emp, '-b', 'linewidth', lw)
plot(vds_emp, coss_emp, '-r', 'linewidth', lw)
plot(vds_emp, crss_emp, '-k', 'linewidth', lw)
hold off
set(gca, 'yscale', 'log')
ylim([1 1e4])
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'fontweight', 'bold', 'fontsize', fsz)

c = zeros(3, n_pts);
for i = 1:n_pts
    vds_cur = vds_sweepvec(i);
    update_netlist(circuit, 'VDS', vds_cur)
    
%     rawfile = [source, '.w'];
%     delete(rawfile)
    
    command = ['"',xyce_path,'" "',circuit,'"'];
%     jsystem(command);
    system(command);
    
    % Read simulation results
    filename = [circuit, '.FD.prn'];
    S = fileread(filename);
    NL = regexp(S, '[\r]');
    headings_raw = S(1:NL(1)-1);
    headings_trm = strtrim(headings_raw);
    headings = strsplit(headings_trm);
    for row = 1:length(NL)-2
        cur_row_raw = S(NL(row):NL(row+1)-1);
        cur_row_trm = strtrim(cur_row_raw); 
        cur_row_cell = strsplit(cur_row_trm); 
    end
    body_raw = S(NL(1)+2:NL(end-1)-1);
    split_data = split(body_raw);
    n = numel(headings); 
    r = numel(split_data)/n; 
    data_num = sprintf('%s ', split_data{:});
    data_col = sscanf(data_num, '%f');
    data_mat = reshape(data_col,n,r);
    data = data_mat';

    vir_idx = strcmpi(headings, 'Re(V(NI))');
    vii_idx = strcmpi(headings, 'Im(V(NI))');
    vor_idx = strcmpi(headings, 'Re(V(NO))');
    voi_idx = strcmpi(headings, 'Im(V(NO))');
    vrr_idx = strcmpi(headings, 'Re(V(NR))');
    vri_idx = strcmpi(headings, 'Im(V(NR))');
    iir_idx = strcmpi(headings, 'Re(I(VISS))');
    iii_idx = strcmpi(headings, 'Im(I(VISS))');
    ior_idx = strcmpi(headings, 'Re(I(VOSS))');
    ioi_idx = strcmpi(headings, 'Im(I(VOSS))');
    irr_idx = strcmpi(headings, 'Re(I(VRSS))');
    iri_idx = strcmpi(headings, 'Im(I(VRSS))');

    vir = data(:,vir_idx);
    vii = data(:,vii_idx);
    vor = data(:,vor_idx);
    voi = data(:,voi_idx);
    vrr = data(:,vrr_idx);
    vri = data(:,vri_idx);
    iir = data(:,iir_idx);
    iii = data(:,iii_idx);
    ior = data(:,ior_idx);
    ioi = data(:,ioi_idx);
    irr = data(:,irr_idx);
    iri = data(:,iri_idx);

    f_idx = strcmpi(headings, 'FREQ');
    f_vec = data(:,f_idx);
    f = unique(f_vec);
    omega = 2*pi*f; 

    Ciss = 1e12*abs((iir.^2+iii.^2)/(vir.*iii-vii.*iir)*1/omega);
    Coss = 1e12*abs((ior.^2+ioi.^2)/(vor.*ioi-voi.*ior)*1/omega);
    Crss = 1e12*abs((irr.^2+iri.^2)/(vrr.*iri-vri.*irr)*1/omega);

    c(1,i) = Ciss;
    c(2,i) = Coss;
    c(3,i) = Crss;
    
    figure(1)
    hold on
    scatter(vds_cur, Ciss, '.b')
    scatter(vds_cur, Coss, '.r')
    scatter(vds_cur, Crss, '.k')
    drawnow
    hold off
    
%     figure(2)
%     hold on
%     scatter(vds_cur, Ciss, '.b')
%     scatter(vds_cur, Coss, '.r')
% %     scatter(vds_cur, Crss, '.k')
%     drawnow
%     hold off
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% % Perform Xyce simulation from specified netlist
% command = ['"',xyce_path,'" ',circuit];
% jsystem(command);
% 
% % Read simulation results
% filename = [circuit, '.FD.prn'];
% S = fileread(filename);
% NL = regexp(S, '[\r]');
% headings_raw = S(1:NL(1)-1);
% headings_trm = strtrim(headings_raw);
% headings = strsplit(headings_trm);
% for row = 1:length(NL)-2
%     cur_row_raw = S(NL(row):NL(row+1)-1);
%     cur_row_trm = strtrim(cur_row_raw); 
%     cur_row_cell = strsplit(cur_row_trm); 
% end
% body_raw = S(NL(1)+2:NL(end-1)-1);
% split_data = split(body_raw);
% n = numel(headings); 
% r = numel(split_data)/n; 
% data_num = sprintf('%s ', split_data{:});
% data_col = sscanf(data_num, '%f');
% data_mat = reshape(data_col,n,r);
% data = data_mat';
% 
% vir_idx = strcmpi(headings, 'Re(V(NI))');
% vii_idx = strcmpi(headings, 'Im(V(NI))');
% vor_idx = strcmpi(headings, 'Re(V(NO))');
% voi_idx = strcmpi(headings, 'Im(V(NO))');
% vrr_idx = strcmpi(headings, 'Re(V(NR))');
% vri_idx = strcmpi(headings, 'Im(V(NR))');
% iir_idx = strcmpi(headings, 'Re(I(VISS))');
% iii_idx = strcmpi(headings, 'Im(I(VISS))');
% ior_idx = strcmpi(headings, 'Re(I(VOSS))');
% ioi_idx = strcmpi(headings, 'Im(I(VOSS))');
% irr_idx = strcmpi(headings, 'Re(I(VRSS))');
% iri_idx = strcmpi(headings, 'Im(I(VRSS))');
% 
% vir = data(:,vir_idx);
% vii = data(:,vii_idx);
% vor = data(:,vor_idx);
% voi = data(:,voi_idx);
% vrr = data(:,vrr_idx);
% vri = data(:,vri_idx);
% iir = data(:,iir_idx);
% iii = data(:,iii_idx);
% ior = data(:,ior_idx);
% ioi = data(:,ioi_idx);
% irr = data(:,irr_idx);
% iri = data(:,iri_idx);
% 
% f_idx = strcmpi(headings, 'FREQ');
% f_vec = data(:,f_idx);
% f = unique(f_vec);
% omega = 2*pi*f; 
% 
% ciss = 1e12*abs((iir.^2+iii.^2)/(vir.*iii-vii.*iir)*1/omega);
% coss = 1e12*abs((ior.^2+ioi.^2)/(vor.*ioi-voi.*ior)*1/omega);
% crss = 1e12*abs((irr.^2+iri.^2)/(vrr.*iri-vri.*irr)*1/omega);
% 
% n = length(vgs); 
% markers = zeros(1,n);
% for i = 1:n
%     cur_vgs = vgs(i);
%     sig_start = find(cur_vgs==vg, 1);
%     markers(i) = sig_start;
% end 
% start_idxs = markers;
% end_idxs = [markers(2:end)-1, length(vd)];
% sig_len = end_idxs(1)-start_idxs(1)+1;
% ids = zeros(n, sig_len);
% for k = 1:n
%     i_idx = start_idxs(k); 
%     f_idx = end_idxs(k); 
%     cur_id = id(i_idx:f_idx);
%     ids(k,:) = cur_id; 
%     if k == 1
%         vds = vd(i_idx:f_idx);
%     end
% end
% 
% 
% lw = 1.2; 
% fsz = 10; 
% figure
% hold on
% % plot(vds_emp, id_emp, 'ok')
% plot(NaN, NaN, '-k', 'linewidth', lw)
% plot(NaN, NaN, ':k', 'linewidth', lw)
% plot(vds_emp, id_emp, '-k', 'linewidth', lw)
% plot(vds, ids, ':k', 'linewidth', lw)
% % % % h1=plot(vds_emp, id_emp, '-k', 'linewidth', 1.2);
% % % % h2=plot(vds, ids, ':k', 'linewidth', 1.2);
% % % % set(get(get(h1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % % % set(get(get(h2(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % plot(vds, ids, '-k', 'linewidth', 1.2)
% % plot(vd, id, '-k', 'linewidth', 1.2)
% % plot(t, v3, '-b', 'linewidth', 1.2)
% % plot(t, v4, '-k', 'linewidth', 1.2)
% hold off
% grid on
% box on
% xlabel('V_D_S (V)')
% ylabel('I_D (A)')
% xlim([0 11])
% ylim([0 50])
% legend({'Datasheet', 'Xyce Simulation'},'location','northwest')
% set(gca, 'fontweight', 'bold', 'fontsize', fsz)
% % legend('Datasheet','Xyce Simulation')
% 
% % legend({'V(2)', 'V(3)', 'V(4)'}) 
% 
% % figure
% % plot(vds, ids, '.k', 'linewidth', lw)


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