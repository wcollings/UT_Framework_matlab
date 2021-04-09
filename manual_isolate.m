% rawdata = read_data('Test4.csv');
rawdata = read_data('Test7.csv');

headings = rawdata(1,:); 
t_loc = strcmpi(headings, 'time');
vds_loc = strcmpi(headings, 'vds');
vgs_loc = strcmpi(headings, 'vgs');
id_loc = strcmpi(headings, 'id');

t = rawdata{2,t_loc};
vgs = rawdata{2,vgs_loc};
vds = rawdata{2,vds_loc};
id = rawdata{2,id_loc};

figure
subplot(3, 1, 1)
plot(t, vds)
subplot(3, 1, 2)
plot(t, id)
subplot(3, 1, 3)
plot(t, vgs)

xmin2 = 2.11e-5;
xmax2 = 2.14e-5;
idx2 = t>=xmin2&t<=xmax2; 
t2 = t(idx2);
vds2 = vds(idx2);
id2 = id(idx2);
vgs2 = vgs(idx2);

figure
subplot(3, 1, 1)
plot(t2, vds2)
axis tight
subplot(3, 1, 2)
plot(t2, id2)
axis tight
subplot(3, 1, 3)
plot(t2, vgs2)
axis tight



% clearvars -except t2 vds2 vgs2 id2
% t_emp = t2; 
% vgs_emp = vgs2; 
% vds_emp = vds2;
% id_emp = id2;
% clearvars -except t_emp id_emp vds_emp vgs_emp



xmin3 = 2.6355e-5;
xmax3 = 2.652e-5;
idx3 = t>=xmin3&t<=xmax3; 
t3 = t(idx3);
vds3 = vds(idx3);
id3 = id(idx3);
vgs3 = vgs(idx3);

figure
subplot(3, 1, 1)
plot(t3, vds3)
axis tight
subplot(3, 1, 2)
plot(t3, id3)
axis tight
subplot(3, 1, 3)
plot(t3, vgs3)
axis tight



% clearvars -except t3 vds3 vgs3 id3
% t_emp = t3; 
% vgs_emp = vgs3; 
% vds_emp = vds3;
% id_emp = id3;
% clearvars -except t_emp id_emp vds_emp vgs_emp