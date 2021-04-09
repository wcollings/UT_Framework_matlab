%%% Read in and format data
load('C:\Users\andre\Documents\Research\Automated SPICE Modeling\MSA12080A_iv_extrap.mat')
% load('C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\Fitting_Data\Static_Data\id_vds_fromdatasheet_25C.mat');
s = 21;
maxs = [nan nan nan nan nan nan 14.8 9 7.2 6.2 6];
for k = 1:size(id_iv,1)
    id_cur = smooth(id_iv(k,:),s);
    id_cur(vds_iv>=maxs(k)) = NaN; 
    id_iv(k,:) = id_cur; 
end
vgs = vgs_iv';
% id_iv = id;
% vds_iv = vds;
x1 = vds_iv'; 
x2 = vds_iv'; 
x3 = vds_iv'; 
x4 = vds_iv'; 
x5 = vds_iv'; 
x6 = vds_iv'; 
x7 = vds_iv'; 
x8 = vds_iv'; 
x9 = vds_iv'; 
x10 = vds_iv'; 
x11 = vds_iv';
y1 = id_iv(1,:)';
y2 = id_iv(2,:)';
y3 = id_iv(3,:)';
y4 = id_iv(4,:)';
y5 = id_iv(5,:)';
y6 = id_iv(6,:)';
y7 = id_iv(7,:)';
y8 = id_iv(8,:)';
y9 = id_iv(9,:)';
y10 = id_iv(10,:)';
y11 = id_iv(11,:)';
id1 = isnan(y1);
id2 = isnan(y2);
id3 = isnan(y3);
id4 = isnan(y4);
id5 = isnan(y5);
id6 = isnan(y6);
id7 = isnan(y7);
id8 = isnan(y8);
id9 = isnan(y9);
id10 = isnan(y10);
id11 = isnan(y11);
x1(id1) = [];
y1(id1) = [];
x2(id2) = [];
y2(id2) = [];
x3(id3) = [];
y3(id3) = [];
x4(id4) = [];
y4(id4) = [];
x5(id5) = [];
y5(id5) = [];
x6(id6) = [];
y6(id6) = [];
x7(id7) = [];
y7(id7) = [];
x8(id8) = [];
y8(id8) = [];
x9(id9) = [];
y9(id9) = [];
x10(id10) = [];
y10(id10) = [];
x11(id11) = [];
y11(id11) = [];
n_sets = length(vgs);
data = cell(n_sets,2);
data{1,1} = x1;
data{2,1} = x2; 
data{3,1} = x3; 
data{4,1} = x4; 
data{5,1} = x5; 
data{6,1} = x6; 
data{7,1} = x7; 
data{8,1} = x8; 
data{9,1} = x9; 
data{10,1} = x10;
data{11,1} = x11; 
data{1,2} = y1;
data{2,2} = y2; 
data{3,2} = y3; 
data{4,2} = y4; 
data{5,2} = y5; 
data{6,2} = y6; 
data{7,2} = y7; 
data{8,2} = y8; 
data{9,2} = y9; 
data{10,2} = y10;
data{11,2} = y11; 

% load C:\Users\andre\Documents\MATLAB\IV_data\Cree
% vgs = [7,9,11,13,15];
% x1 = Cree_vgs7(:,1);
% x2 = Cree_vgs9(:,1);
% x3 = Cree_vgs11(:,1);
% x4 = Cree_vgs13(:,1);
% x5 = Cree_vgs15(:,1);
% y1 = Cree_vgs7(:,2);
% y2 = Cree_vgs9(:,2);
% y3 = Cree_vgs11(:,2);
% y4 = Cree_vgs13(:,2);
% y5 = Cree_vgs15(:,2);
% n_sets = length(vgs);
% data = cell(n_sets,2);
% data{1,1} = x1;
% data{2,1} = x2;
% data{3,1} = x3;
% data{4,1} = x4;
% data{5,1} = x5;
% data{1,2} = y1;
% data{2,2} = y2;
% data{3,2} = y3;
% data{4,2} = y4;
% data{5,2} = y5;

lw=1.2;
figure(1)
hold on
plot(x1,y1,'-k','linewidth',lw)
plot(x2,y2,'-k','linewidth',lw)
plot(x3,y3,'-k','linewidth',lw)
plot(x4,y4,'-k','linewidth',lw)
plot(x5,y5,'-k','linewidth',lw)
plot(x6,y6,'-k','linewidth',lw)
plot(x7,y7,'-k','linewidth',lw)
plot(x8,y8,'-k','linewidth',lw)
plot(x9,y9,'-k','linewidth',lw)
plot(x10,y10,'-k','linewidth',lw)
plot(x11,y11,'-k','linewidth',lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('I_D (A)')
text(8.5,5.4,'V_G_S = 7 V')
text(9,10.2,'9 V')
text(8.2,15.4,'11 V')
text(7.5,18.4,'13 V')
text(7,19.5,'15 V')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
dx1 = mean([x1(1:end-1),x1(2:end)],2);
dx2 = mean([x2(1:end-1),x2(2:end)],2);
dx3 = mean([x3(1:end-1),x3(2:end)],2);
dx4 = mean([x4(1:end-1),x4(2:end)],2);
dx5 = mean([x5(1:end-1),x5(2:end)],2);
dx6 = mean([x6(1:end-1),x6(2:end)],2);
dx7 = mean([x7(1:end-1),x7(2:end)],2);
dx8 = mean([x8(1:end-1),x8(2:end)],2);
dx9 = mean([x9(1:end-1),x9(2:end)],2);
dx10 = mean([x10(1:end-1),x10(2:end)],2);
dx11 = mean([x11(1:end-1),x11(2:end)],2);
dy1 = diff(y1)./diff(x1); 
dy2 = diff(y2)./diff(x2); 
dy3 = diff(y3)./diff(x3); 
dy4 = diff(y4)./diff(x4); 
dy5 = diff(y5)./diff(x5); 
dy6 = diff(y6)./diff(x6);
dy7 = diff(y7)./diff(x7);
dy8 = diff(y8)./diff(x8);
dy9 = diff(y9)./diff(x9);
dy10 = diff(y10)./diff(x10);
dy11 = diff(y11)./diff(x11);
% n_pts = 100; 
% x1i = linspace(x1(1),x1(end),n_pts);
% x2i = linspace(x2(1),x2(end),n_pts);
% x3i = linspace(x3(1),x3(end),n_pts);
% x4i = linspace(x4(1),x4(end),n_pts);
% x5i = linspace(x5(1),x5(end),n_pts);
% y1i = pchip(x1, y1, x1i);
% y2i = pchip(x2, y2, x2i);
% y3i = pchip(x3, y3, x3i);
% y4i = pchip(x4, y4, x4i);
% y5i = pchip(x5, y5, x5i);
% dx1 = mean([x1i(1:end-1);x1i(2:end)],1);
% dx2 = mean([x2i(1:end-1);x2i(2:end)],1);
% dx3 = mean([x3i(1:end-1);x3i(2:end)],1);
% dx4 = mean([x4i(1:end-1);x4i(2:end)],1);
% dx5 = mean([x5i(1:end-1);x5i(2:end)],1);
% dy1 = diff(y1i)./diff(x1i); 
% dy2 = diff(y2i)./diff(x2i); 
% dy3 = diff(y3i)./diff(x3i); 
% dy4 = diff(y4i)./diff(x4i); 
% dy5 = diff(y5i)./diff(x5i); 
figure(2)
hold on
plot(dx1,dy1,'-k','linewidth',lw)
plot(dx2,dy2,'-k','linewidth',lw)
plot(dx3,dy3,'-k','linewidth',lw)
plot(dx4,dy4,'-k','linewidth',lw)
plot(dx5,dy5,'-k','linewidth',lw)
plot(dx6,dy6,'-k','linewidth',lw)
plot(dx7,dy7,'-k','linewidth',lw)
plot(dx8,dy8,'-k','linewidth',lw)
plot(dx9,dy9,'-k','linewidth',lw)
plot(dx10,dy10,'-k','linewidth',lw)
plot(dx11,dy11,'-k','linewidth',lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('Conductance (A/V)')
% text(8.5,0.375,'V_G_S = 7 V')
% text(8.1,0.75,'9 V')
% text(7.35,1.4,'11 V')
% text(6.85,2,'13 V')
% text(6.3,2.3,'15 V')
set(gca, 'fontweight', 'bold', 'fontsize', 10)

%%% Generate forward conduction model
fit_coeff = zeros(n_sets,3);
figure(3)
figure(4)
alpha_max = 0; 
beta_max = inf;
zeta_min = 0; 
for dataset = 1:n_sets
    x = data{dataset,1};
    y = data{dataset,2};
    dy = diff(y)./diff(x);
    dx = mean([x(1:end-1),x(2:end)],2);
    dfun = @(x,xdata)x(1)*(atan(x(2)*xdata)-pi/2)+x(3);
    alpha = -0.5; 
    beta = 0.5; 
    zeta = 0.5; 
    x0 = [alpha, beta, zeta];
%     x0 = [-0.5, 0.5, 0.5];
%     lb = [-inf, 0, 0];
%     ub = [0, inf, inf];
    lb = [-inf, 0, zeta_min];
    ub = [alpha_max, beta_max, inf];
    dfit = lsqcurvefit(dfun,x0,dx,dy,lb,ub);
    dy_fit = dfun(dfit,dx);
    alpha_max = dfit(1);
    beta_min = dfit(2);
    zeta_min = dfit(3);
%     figure
    figure(3)
    hold on
%     plot(x_new,dy)
    scatter(dx,dy,'.k')
    plot(dx,dy_fit,':k','linewidth',lw)
    hold off
    fit_coeff(dataset,:) = dfit; 
    %%% idfun = -(a*log(b^2*x^2+1)/(2*b)+a*x*(atan(b*x)-pi/2)+z*x
    idfun = @(x,xdata)-(x(1).*log(x(2)^2.*xdata.^2+1))/(2.*x(2))+x(1).*xdata.*(atan(x(2).*xdata)-pi/2)+x(3).*xdata;
    id = idfun(dfit,x);
    x_ext = linspace(0, 100, 1000);
    id_ext = idfun(dfit, x_ext);
%     figure
    figure(4)
    hold on
    scatter(x,y,'.k')
    plot(x,id,'-k','linewidth',lw)
%     plot(x_ext,id_ext)
    hold off
end
figure (3)
grid on
box on
xlabel('V_D_S (V)')
ylabel('Conductance (A/V)')
% text(8.5,0.375,'V_G_S = 7 V')
% text(8.1,0.75,'9 V')
% text(7.35,1.4,'11 V')
% text(6.85,2,'13 V')
% text(6.3,2.3,'15 V')
set(gca, 'fontweight', 'bold', 'fontsize', 10)

figure(4)
grid on
box on
xlabel('V_D_S (V)')
ylabel('I_D (A)')
% text(8.5,5.4,'V_G_S = 7 V')
% text(9,10.2,'9 V')
% text(8.2,15.4,'11 V')
% text(7.5,18.4,'13 V')
% text(7,19.5,'15 V')
set(gca, 'fontweight', 'bold', 'fontsize', 10)

% figure
% subplot(3, 1, 1)
% scatter(vgs,fit_coeff(:,1))
% xlabel('Vgs (V)')
% ylabel('c1')
% subplot(3, 1, 2)
% scatter(vgs,fit_coeff(:,2))
% xlabel('Vgs (V)')
% ylabel('c2')
% subplot(3, 1, 3)
% scatter(vgs,fit_coeff(:,3))
% xlabel('Vgs (V)')
% ylabel('c3')

center = fit_coeff(2,:); 
sweep_perc = 50;
n_vals = 5;
% x = linspace(0,100,1000);
titles = {'\alpha', '\beta', '\zeta'};
% figure('units', 'normalized', 'outerposition', [0, 0, 0.8, 0.6])
sweep_start = [13, 11.55, 7.25];
sweep_end = [3.75, 6.25, 9.6];
for i = 1:length(center)
    coeff = center(i);
    sweep_i = coeff - abs(sweep_perc/100*coeff);
    sweep_f = coeff + abs(sweep_perc/100*coeff);
    sweepvec = linspace(sweep_i, sweep_f, n_vals);
    node = center;
%     subplot(1, 3, i)
    figure(4+i)
    hold on
    for k = 1:n_vals
        node(i) = sweepvec(k);
        curve = idfun(node,x);
        if k == (n_vals + 1)/2
            plot(x, curve, '-k', 'linewidth', lw)
        else
            plot(x, curve, ':k', 'linewidth', lw)
        end
    end
    text(6, sweep_start(i), '-50%')
    text(6, sweep_end(i), '+50%')
    hold off
    cur_title = titles{i};
    title(cur_title, 'fontweight', 'bold', 'fontsize', 14)
    xlim([x(1), x(end)])
    ylim([0 18])
    xlabel('V_D_S (V)')
    ylabel('I_D (A)')
    set(gca, 'fontweight', 'bold', 'fontsize', 10)
    grid on
    box on
%     ylim([-100 200])
end

figure(8)
alpha_vals = fit_coeff(:,1);
beta_vals = fit_coeff(:,2);
zeta_vals = fit_coeff(:,3);
subplot(3, 1, 1)
scatter(vgs, alpha_vals, 'k', 'filled')
ylabel('\alpha')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on
subplot(3, 1, 2)
scatter(vgs, beta_vals, 'k', 'filled')
ylabel('\beta')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on
subplot(3, 1, 3)
scatter(vgs, zeta_vals, 'k', 'filled')
ylabel('\zeta')
xlabel('V_G_S (V)')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on

% close all
sigfun = @(x,xdata)x(1)./(1+exp(-x(2)*xdata+x(3)))+x(4);
x = 0:0.1:20;

z0 = [0.4, 1, 9, 0];
zlb = [0, 0, -inf, 0];
zub = [inf, inf, inf, inf];
zfit = lsqcurvefit(sigfun, z0, vgs, zeta_vals', zlb, zub);
z = sigfun(zfit, x);
% figure
% hold on
% plot(x,z)
% scatter(vgs, zeta_vals)
% hold off
a0 = [0.4, 1, 9, 0];
alb = [-inf, 0, -inf, -inf];
aub = [0, inf, inf, 0];
afit = lsqcurvefit(sigfun, a0, vgs, alpha_vals', alb, aub);
a = sigfun(afit, x);
% figure
% hold on
% plot(x,a)
% scatter(vgs, alpha_vals)
% hold off
bsigfun = @(x,xdata)-x(1)./(1+exp(-x(2)*xdata+x(3)))+x(1)+x(4);
b0 = [0.4, 1, 9, 0];
blb = [0, 0, -inf, 0];
bub = [inf, inf, inf, inf];
bfit = lsqcurvefit(bsigfun, b0, vgs, beta_vals', blb, bub);
b = bsigfun(bfit, x);
% figure
% hold on
% plot(x,b)
% scatter(vgs, beta_vals)
% hold off
figure(9)
alpha_vals = fit_coeff(:,1);
beta_vals = fit_coeff(:,2);
zeta_vals = fit_coeff(:,3);
subplot(3, 1, 1)
hold on
scatter(vgs, alpha_vals, 'k', 'filled')
plot(x,a,'-k','linewidth',lw)
hold off
ylabel('\alpha')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on
xlim([7,15])
subplot(3, 1, 2)
hold on
scatter(vgs, beta_vals, 'k', 'filled')
plot(x,b,'-k','linewidth',lw)
hold off
ylabel('\beta')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on
xlim([7,15])
subplot(3, 1, 3)
hold on
scatter(vgs, zeta_vals, 'k', 'filled')
plot(x,z,'-k','linewidth',lw)
hold off
ylabel('\zeta')
xlabel('V_G_S (V)')
set(gca, 'fontweight', 'bold', 'fontsize', 10)
grid on
box on
xlim([7,15])


%%% Plot intial fit using individually fit coefficients
vds = 0:0.1:10;
if isrow(vgs)
    vgs = vgs';
end
alpha = sigfun(afit, vgs);
beta = bsigfun(bfit, vgs);
zeta = sigfun(zfit, vgs);
idfitfun = @(x,xdata)-(alpha.*log(beta.^2.*xdata.^2+1))./(2.*beta)+alpha.*xdata.*(atan(beta.*xdata)-pi/2)+zeta.*xdata;
id_fit_individual = idfitfun([],vds);
figure
hold on
plot(vds, id_fit_individual,'-k')
scatter(x1,y1,'k')
scatter(x2,y2,'k')
scatter(x3,y3,'k')
scatter(x4,y4,'k')
scatter(x5,y5,'k')

vgs_ext = 0:1:20;
vgs_ext = vgs_ext';
alpha = sigfun(afit, vgs_ext);
beta = bsigfun(bfit, vgs_ext);
zeta = sigfun(zfit, vgs_ext);
c0 = [alpha, beta, zeta];
idfitfun = @(x,xdata)-(alpha.*log(beta.^2.*xdata.^2+1))./(2.*beta)+alpha.*xdata.*(atan(beta.*xdata)-pi/2)+zeta.*xdata;
id_fit_individual = idfitfun([],vds);
vgs_data_idx = find(ismember(vgs_ext,[vgs]));
figure
hold on
plot(vds, id_fit_individual,':k','linewidth',lw)
plot(vds, id_fit_individual(vgs_data_idx,:),'-k','linewidth',lw)
scatter(x1,y1,'k')
scatter(x2,y2,'k')
scatter(x3,y3,'k')
scatter(x4,y4,'k')
scatter(x5,y5,'k')
hold off
xlabel('V_D_S (V)')
ylabel('I_D (A)')
grid on
box on
set(gca,'fontweight','bold','fontsize',10)
ylim([0 20])

%%% Perform simulataneous fitting of all of the fit parameters. 
k0 = [afit, bfit, zfit];
flb = [alb, blb, zlb];
fub = [aub, bub, zub];
% alpha = (k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4))
% beta = (k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8))
% zeta = (k(9)./(1+exp(-k(10)*vgs+k(11)))+k(12))
% full_fun=@(k,xdata)-((k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4)).*log((k(5).*(1-1/(1+exp(-k(6)*vgs+k(7))))+k(8)).^2.*xdata.^2+1))./(2.*(k(5).*(1-1/(1+exp(-k(6)*vgs+k(7))))+k(8)))+(k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4)).*xdata.*(atan((k(5).*(1-1/(1+exp(-k(6)*vgs+k(7))))+k(8)).*xdata)-pi/2)+(k(9)./(1+exp(-k(10)*vgs+k(11)))+k(12)).*xdata;
% full_fun=@(k,xdata)-((k(1)./(1+exp(-k(2).*vgs+k(3)))+k(4)).*log((k(5).*(1-1./(1+exp(-k(6).*vgs+k(7))))+k(8)).^2*xdata.^2+1))./(2.*(k(5).*(1-1/(1+exp(-k(6).*vgs+k(7))))+k(8)))+(k(1)./(1+exp(-k(2).*vgs+k(3)))+k(4))*xdata*(atan((k(5).*(1-1/(1+exp(-k(6).*vgs+k(7))))+k(8))*xdata)-pi/2)+(k(9)./(1+exp(-k(10).*vgs+k(11)))+k(12))*xdata;

% m1 = (-(alpha.*log(beta.^2*xdata.^2+1))./(2*beta))
% m2 = (alpha*xdata.*(atan(beta*xdata)-pi/2))
% m3 = (zeta*xdata);
% m1 = (-((k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4)).*log((k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8)).^2*xdata.^2+1))./(2*(k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8))));
% m2 = ((k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4))*xdata.*(atan((k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8))*xdata)-pi/2));
% m3 = ((k(9)./(1+exp(-k(10)*vgs+k(11)))+k(12))*xdata);
% m4 = m1+m2+m3;

% fun = @(k,xdata)(-((k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4)).*log((k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8)).^2*xdata.^2+1))./(2*(k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8))))+((k(1)./(1+exp(-k(2)*vgs+k(3)))+k(4))*xdata.*(atan((k(5)*(1-1./(1+exp(-k(6)*vgs+k(7))))+k(8))*xdata)-pi/2))+((k(9)./(1+exp(-k(10)*vgs+k(11)))+k(12))*xdata);
% initial_sol = full_fun(k0,vds);
