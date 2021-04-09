%%% Set Up Workspace
% Add nested directories to path
folder = fileparts(which(mfilename));
addpath(genpath(folder));

warning('off','all')

rawdataCV = readtable('EPC2022_caps.csv');
% v_cv = rawdataCV.VDS;
% ciss = rawdataCV.Ciss;
% coss = rawdataCV.Coss;
% crss = rawdataCV.Crss;
ds=rawdataCV;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % file = 'cv_testbed';
file = 'cv_testbed';
% % file = 'cv_testbed_wolf';


vds_start = 0;
vds_end = log10(100);
n_pts = 30;
% n_pts = 1000;

% vds_start = 0;
% vds_end = 50;
% % n_pts = 30;
% n_pts = 2000;

% v_cv = 5; 

lw = 1.2;
fsz=10;
figure(1)
clf
hold on
plot(ds.viss, ds.ciss, '-b', 'linewidth', lw)
plot(ds.voss, ds.coss, '-r', 'linewidth', lw)
plot(ds.vrss, ds.crss, '-k', 'linewidth', lw)
hold off
% set(gca, 'xscale', 'log'
set(gca, 'yscale', 'log')
% ylim([1e1 1e4])
ylim([1 1e4])
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'fontweight', 'bold', 'fontsize', fsz)
curr=pwd;
path = append(curr, "\SPICE_Files\");
source = append(path,file);
ltspicepath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe';
circuit = append(source, '.asc');

vds_sweepvec = logspace(vds_start, vds_end, n_pts);
% vds_start = -5;
% vds_end = 5;
% n_pts = 100;
% vds_sweepvec = linspace(vds_start, vds_end, n_pts);

f = 1e6;
w = 2*pi*f;

c = zeros(3, n_pts);
% figure
for i = 1:n_pts
    vds_cur = vds_sweepvec(i);
    setvds(circuit, vds_cur)
    
%     rawfile = [source, '.w'];
%     delete(rawfile)
    rawfile = append(source, '.raw');
    %delete(rawfile);    

    command = ['"',ltspicepath,'"', ' -run -b',' ',circuit];
    jsystem(strjoin(command, ''));
    
    rawdata = LTspice2Matlab("SPICE_Files\cv_testbed.raw");
    
%     rawfile = [source, '.raw'];
    
    headings = rawdata.variable_name_list;    
    Iiss_idx = strcmpi(headings, 'I(V_iss)');
    Ioss_idx = strcmpi(headings, 'I(V_oss)');
    Irss_idx = strcmpi(headings, 'I(V_rss)');
    Viss_idx = strcmpi(headings, 'V(stim1)');
    Voss_idx = strcmpi(headings, 'V(stim2)');
    Vrss_idx = strcmpi(headings, 'V(stim3)');
    
    data = rawdata.variable_mat;
    Iiss = data(Iiss_idx);
    Ioss = data(Ioss_idx);
    Irss = data(Irss_idx);
    Viss = data(Viss_idx);
    Voss = data(Voss_idx);
    Vrss = data(Vrss_idx);
    
    Ciss = 1e12*abs((real(Iiss)^2+imag(Iiss)^2)/(real(Viss)*imag(Iiss)-imag(Viss)*real(Iiss))*1/w);
    Coss = 1e12*abs((real(Ioss)^2+imag(Ioss)^2)/(real(Voss)*imag(Ioss)-imag(Voss)*real(Ioss))*1/w);
    Crss = 1e12*abs((real(Irss)^2+imag(Irss)^2)/(real(Vrss)*imag(Irss)-imag(Vrss)*real(Irss))*1/w);
%     Ciss = abs((real(Iiss)^2+imag(Iiss)^2)/(real(Viss)*imag(Iiss)-imag(Viss)*real(Iiss))*1/w);
%     Coss = abs((real(Ioss)^2+imag(Ioss)^2)/(real(Voss)*imag(Ioss)-imag(Voss)*real(Ioss))*1/w);
%     Crss = abs((real(Irss)^2+imag(Irss)^2)/(real(Vrss)*imag(Irss)-imag(Vrss)*real(Irss))*1/w);
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
    
end

% figure
% plot(vds_sweepvec, c)
% set(gca, 'yscale', 'log')
% set(gca, 'xscale', 'log')
% legend({'C_I_S_S','C_O_S_S','C_R_S_S'})
% ylim([10 1e4])

% ciss_sim = c(1,:);
% coss_sim = c(2,:);
% crss_sim = c(3,:);
% cdg_sim = crss_sim;
% cds_sim = coss_sim - crss_sim;
% cgs_sim = ciss_sim - crss_sim;
% figure
% hold on
% plot(vds_sweepvec, cgs_sim)
% plot(vds_sweepvec, cds_sim)
% plot(vds_sweepvec, cdg_sim)
% hold off
% set(gca, 'yscale', 'log')
% legend({'C_G_S','C_D_S','C_D_G'})
% set(gca, 'xscale', 'log')
% ylim([10 1e4])


% % c_const = 1e-12;
% clearvars -except vds_sweepvec cgs_sim cds_sim cdg_sim
% % plot(vds_sweepvec, cgs_sim)
% % plot(vds_sweepvec, cds_sim)
% % plot(vds_sweepvec, cdg_sim)


function setvds(netlist, vds)
    fcell = {};
    fid = fopen(netlist,'r');
    nextline = fgetl(fid);
    n_line = 1;
    while nextline ~= -1
        targidx = regexpi(nextline, '.param VDS');
        if ~isempty(targidx)
            curline = strsplit(nextline);
            vds_loc = strcmpi(curline, 'VDS');
            targ_id = find(vds_loc) + 2;
            curline{targ_id} = num2str(vds);
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