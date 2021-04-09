
%%%%%% Read in data and set preferences
% cvsource = 'C:\Users\andre\Documents\Research\SiC_Modeling\CREE Model\cvdata.mat';
% load(cvsource)


cvsource = 'C:\Users\andre\Documents\Research\Automated SPICE Modeling\UT_Framework_3.0\Fitting_Data\Static_Data\NewData_CVData_pf.csv';
rawcv = read_data(cvsource);
sigs = {'C_I_S_S', 'C_O_S_S', 'C_R_S_S'};
vds = rawcv{2,1}; 
c = [rawcv{2,2};rawcv{2,3};rawcv{2,4}];

lw = 1.2;
fsz = 10; 

vds_new = logspace(-1, log10(vds(end)), 1000);
% c = interp1(vds, c, vds_new, 'linear', 'extrap');
c = pchip(vds, c, vds_new);
vds = vds_new; 


%%%%%% Plot data
figure(1)
plot(vds, c, 'linewidth', lw)
legend(sigs)
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1 1e4])

crss_idx = strcmpi(sigs, 'C_R_S_S');
ciss_idx = strcmpi(sigs, 'C_I_S_S');
coss_idx = strcmpi(sigs, 'C_O_S_S');
crss = c(crss_idx,:);
ciss = c(ciss_idx,:);
coss = c(coss_idx,:);
cgd = crss;
cgs = ciss - cgd;
cds = coss - cgd; 

figure(2)
hold on
plot(vds, cgs, '-r', 'linewidth', lw)
plot(vds, cds, '-b', 'linewidth', lw)
plot(vds, cgd, '-k', 'linewidth', lw)
hold off
legend({'C_G_S', 'C_D_S', 'C_G_D'})
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1e0 1e4])






%%%%%% Perform "traditional fitting" for comparison
%%% For standard form, C = Co/((1+a*vds)^chi)
cogs = pchip(vds, cgs, 0);
cods = pchip(vds, cds, 0);
cogd = pchip(vds, cgd, 0); 

%%% Direct fitting with exponential form
efun = @(x,xdata)x(1)./((1+x(2)*xdata).^x(3));
egs0 = [cogs, .05, .5];
egspar = lsqcurvefit(efun, egs0, vds, cgs);
cgse = efun(egspar, vds); 

eds0 = [cods, .05, .5];
edspar = lsqcurvefit(efun, eds0, vds, cds);
cdse = efun(edspar, vds); 

egd0 = [cogd, .05, .5];
egdpar = lsqcurvefit(efun, egd0, vds, cgd);
cgde = efun(egdpar, vds); 

figure(3)
hold on
plot(vds, cgs, '-r', 'linewidth', lw)
plot(vds, cds, '-b', 'linewidth', lw)
plot(vds, cgd, '-k', 'linewidth', lw)
plot(vds, cgse, ':r', 'linewidth', lw)
plot(vds, cdse, ':b', 'linewidth', lw)
plot(vds, cgde, ':k', 'linewidth', lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1e0 1e4])
title('Direct Fitting w/ Exponential')

%%% Skipping initial data points
fit_idx = 18:length(vds); 
egspar2 = lsqcurvefit(efun, egs0, vds(fit_idx), cgs(fit_idx));
cgse2 = efun(egspar2, vds); 
edspar2 = lsqcurvefit(efun, eds0, vds(fit_idx), cds(fit_idx));
cdse2 = efun(edspar2, vds); 
egdpar2 = lsqcurvefit(efun, egd0, vds(fit_idx), cgd(fit_idx));
cgde2 = efun(egdpar2, vds); 

figure(4)
hold on
% plot(vds(fit_idx), cgs(fit_idx), '-r', 'linewidth', lw)
% plot(vds(fit_idx), cds(fit_idx), '-b', 'linewidth', lw)
% plot(vds(fit_idx), cgd(fit_idx), '-k', 'linewidth', lw)
% plot(vds(fit_idx), cgse2(fit_idx), ':r', 'linewidth', lw)
% plot(vds(fit_idx), cdse2(fit_idx), ':b', 'linewidth', lw)
% plot(vds(fit_idx), cgde2(fit_idx), ':k', 'linewidth', lw)
plot(vds, cgs, '-r', 'linewidth', lw)
plot(vds, cds, '-b', 'linewidth', lw)
plot(vds, cgd, '-k', 'linewidth', lw)
plot(vds, cgse2, ':r', 'linewidth', lw)
plot(vds, cdse2, ':b', 'linewidth', lw)
plot(vds, cgde2, ':k', 'linewidth', lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1e0 1e4])
title('Exponential Fit Omitting Low V_D_S')




%%% Exponential definition fitting
gamma = 0.5;
alpha = 1; 
xgs = zeros(size(vds)); 
xds = zeros(size(vds)); 
xgd = zeros(size(vds)); 
for k = 1:length(vds)
%     v = vds(k);
    xgs(k) = 1/alpha*((cogs/cgs(k))^(1/gamma)-1);
    xds(k) = 1/alpha*((cods/cds(k))^(1/gamma)-1);
    xgd(k) = 1/alpha*((cogd/cgd(k))^(1/gamma)-1);
end

% figure(3)
% hold on
% plot(vds, xgs)
% plot(vds, xds)
% plot(vds, xgd)
% hold off

xfun = @(x,xdata)x(1)*tanh(x(2)*xdata);
xgs0 = [xgs(end) 0.05];
xgspar = lsqcurvefit(xfun, xgs0, vds, xgs);
xgs_fit = xfun(xgspar, vds);

xds0 = [xds(end) 0.05];
xdspar = lsqcurvefit(xfun, xds0, vds, xds);
xds_fit = xfun(xdspar, vds);

xgd0 = [xgd(end) 0.05];
xgdpar = lsqcurvefit(xfun, xgd0, vds, xgd);
xgd_fit = xfun(xgdpar, vds);

figure(5)
hold on
plot(vds, xgs, '-r', 'linewidth', lw)
plot(vds, xds, '-b', 'linewidth', lw)
plot(vds, xgd, '-k', 'linewidth', lw)
plot(vds, xgs_fit, ':r', 'linewidth', lw)
plot(vds, xds_fit, ':b', 'linewidth', lw)
plot(vds, xgd_fit, ':k', 'linewidth', lw)
hold off

% compfun = @(x,xdata)x(1)./((1+x(2)*vds).^x(3));
% xfun = @(x,xdata)x(1)*tanh(x(2)*xdata);
cgs_comp = cogs./((1+alpha*(xgs_fit)).^gamma);
cds_comp = cods./((1+alpha*(xds_fit)).^gamma);
cgd_comp = cogd./((1+alpha*(xgd_fit)).^gamma);

figure(6)
hold on
plot(vds, cgs, '-r', 'linewidth', lw)
plot(vds, cds, '-b', 'linewidth', lw)
plot(vds, cgd, '-k', 'linewidth', lw)
plot(vds, cgs_comp, ':r', 'linewidth', lw)
plot(vds, cds_comp, ':b', 'linewidth', lw)
plot(vds, cgd_comp, ':k', 'linewidth', lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1e0 1e4])
title('C(f(x)): f(x) = \beta_1\bullettanh(\beta_2\bulletV_D_S)')

% % cgd_comp = cogd./((1+alpha*(6.5801e+03*tanh(0.0028*vds))).^gamma);
% % fullfun = @(x,xdata)cogd./((1+alpha*(6.5801e+03*tanh(0.0028*vds))).^gamma);
% fullfun = @(x,xdata)x(1)./((1+x(2)*(x(3)*tanh(x(4)*xdata))).^x(5));
% x0 = [cogd alpha xgdpar gamma];
% % % fullpar = lsqcurvefit(fullfun, x0, vds, cgd);
% % % fullfit = fullfun(fullpar, vds);
% % figure
% % hold on
% % plot(vds, cgd)
% % plot(vds, fullfun(x0,vds))
% % % plot(vds, fullfit)
% % hold off






%%% Investigate manipulating low-Vds behavior
% xfun = @(x,xdata)x(1)*tanh(x(2)*xdata);
% xgs0 = [xgs(end) 0.05];
% xgspar = lsqcurvefit(xfun, xgs0, vds, xgs);
% xgs_fit = xfun(xgspar, vds);
% 
% xds0 = [xds(end) 0.05];
% xdspar = lsqcurvefit(xfun, xds0, vds, xds);
% xds_fit = xfun(xdspar, vds);
% 
% xgd0 = [xgd(end) 0.05];
% xgdpar = lsqcurvefit(xfun, xgd0, vds, xgd);
% xgd_fit = xfun(xgdpar, vds);
% 
% figure(5)
% hold on
% plot(vds, xgs, '-r', 'linewidth', lw)
% plot(vds, xds, '-b', 'linewidth', lw)
% plot(vds, xgd, '-k', 'linewidth', lw)
% plot(vds, xgs_fit, ':r', 'linewidth', lw)
% plot(vds, xds_fit, ':b', 'linewidth', lw)
% plot(vds, xgd_fit, ':k', 'linewidth', lw)
% hold off

xgs_gain = xgs./xgs_fit; 
xds_gain = xds./xds_fit; 
xgd_gain = xgd./xgd_fit; 

figure(7)
hold on
plot(vds, xgs_gain)
plot(vds, xds_gain)
plot(vds, xgd_gain)
hold off


xsfun = @(x,xdata)x(1)*tanh(x(2)*(xdata-x(3)))+x(4);
xs0 = [.5 1 20 1.1];

kgspar = lsqcurvefit(xsfun,xs0,vds,xgs_gain);
kgs = xsfun(kgspar,vds); 

kdspar = lsqcurvefit(xsfun,xs0,vds,xds_gain);
kds = xsfun(kdspar,vds); 

kgdpar = lsqcurvefit(xsfun,xs0,vds,xgd_gain);
kgd = xsfun(kgdpar,vds); 

figure(8)
hold on
plot(vds, xgs, '-r', 'linewidth', lw)
plot(vds, xds, '-b', 'linewidth', lw)
plot(vds, xgd, '-k', 'linewidth', lw)
plot(vds, xgs_fit.*kgs, ':r', 'linewidth', lw)
plot(vds, xds_fit.*kds, ':b', 'linewidth', lw)
plot(vds, xgd_fit.*kgd, ':k', 'linewidth', lw)
hold off





%%%%%% Apply sigmoid gain factor
fullfun = @(x,xdata)x(1)./((1+x(2)*(x(3)*tanh(x(4)*xdata))).^x(5));
fgs = [cogs alpha xgspar gamma];
cgs_full = fullfun(fgs, vds.*kgs);
fds = [cods alpha xdspar gamma];
cds_full = fullfun(fds, vds.*kds);
fgd = [cogd alpha xgdpar gamma];
cgd_full = fullfun(fgd, vds.*kgd);
figure(9)
hold on
plot(vds, cgs, '-r', 'linewidth', lw)
plot(vds, cds, '-b', 'linewidth', lw)
plot(vds, cgd, '-k', 'linewidth', lw)
plot(vds, cgs_full, ':r', 'linewidth', lw)
plot(vds, cds_full, ':b', 'linewidth', lw)
plot(vds, cgd_full, ':k', 'linewidth', lw)
hold off
grid on
box on
xlabel('V_D_S (V)')
ylabel('Capacitance (pF)')
set(gca, 'YScale', 'log')
set(gca,'fontweight','bold','fontsize',fsz)
ylim([1e0 1e4])










cogd = 1.454e-09;
xgd1 = 15041;
xgd2 = 0.0025;
ggd1 = 0.6192;
ggd2 = 0.4130;
ggd3 = 3.5665;
ggd4 = 0.5415;
agd = 1;
chigd = 0.5;

xgd = xgd1*tanh(xgd2*vds);
cgdmod = cogd*(1+agd+xgd).^(-chigd);


% % %%% Full function including k in one line
% % fullfun = @(x,xdata)x(1)./((1+x(2)*(x(3)*tanh(x(4)*xdata))).^x(5));


% %%%%%% Take additional datapoints to smooth behavior
% vds_new = logspace(log10(vds(1)), log10(vds(end)), 1000);
% kgs_new = xsfun(kgspar,vds_new); 
% kds_new = xsfun(kdspar,vds_new); 
% kgd_new = xsfun(kgdpar,vds_new); 
% cgs_full_new = fullfun(fgs, vds_new.*kgs_new);
% cds_full_new = fullfun(fds, vds_new.*kds_new);
% cgd_full_new = fullfun(fgd, vds_new.*kgd_new);
% figure(10)
% hold on
% plot(vds, cgs, '-r', 'linewidth', lw)
% plot(vds, cds, '-b', 'linewidth', lw)
% plot(vds, cgd, '-k', 'linewidth', lw)
% plot(vds_new, cgs_full_new, ':r', 'linewidth', lw)
% plot(vds_new, cds_full_new, ':b', 'linewidth', lw)
% plot(vds_new, cgd_full_new, ':k', 'linewidth', lw)
% hold off
% grid on
% box on
% xlabel('V_D_S (V)')
% ylabel('Capacitance (pF)')
% set(gca, 'YScale', 'log')
% set(gca,'fontweight','bold','fontsize',fsz)
% ylim([1e0 1e4])







% %%% Implement similar process a second time
% % xgdpar
% % alpha 
% % gamma 
% % cogd
% 
% v = zeros(size(vds)); 
% for k = 1:length(vds)
%     v(k) = 1/xgdpar(2)*atanh(((cogd/cgd(k))*(1/gamma)-1)/(alpha*xgdpar(1)));
% end
% 
% % % % vfun = @(x,xdata)x(1)./(1+x(2)*exp(-x(3)*(xdata-x(4))));
% % % % vgd0 = [v(end) 1 .02 0];
% % % % vgdpar = lsqcurvefit(vfun, vgd0, vds, v);
% % % % vgd_fit = vfun(vgdpar, vds);
% % % 
% % % vgd0 = [v(end) 0.05];
% % % vgdpar = lsqcurvefit(xfun, vgd0, vds, v);
% % % vgd_fit = xfun(vgdpar, vds);
% % % 
% % % figure
% % % hold on
% % % plot(vds, v)
% % % plot(vds, vgd_fit)
% % % hold off
% % % 
% % % % fullfun = @(x,xdata)x(1)./((1+x(2)*(x(3)*tanh(x(4)*vds))).^x(5));
% % % % x0 = [cogd alpha xgdpar gamma];
% % % 
% % % cgd_2 = fullfun(x0,vgd_fit);
% % % 
% % % figure
% % % hold on
% % % plot(vds, cgd)
% % % plot(vds, cgd_2)
% % % hold off
% 
% % %%% hyperbolic tangent fitting
% % % plot(vds, (cogd - cgd(end)).*(tanh(-.05*vds)+1) + cgd(end))
% % alpha_gd = cogd-cgd(end);
% % gamma_gd = cgd(end);
% % beta_gd = 0.05;
% % 
% % xgd = zeros(size(vds)); 
% % for k = 1:length(vds)
% %     v = vds(k);
% %     xgd(k) = -1/beta_gd*atanh((cgd(k)-gamma_gd)/alpha_gd-1);
% % end
% % 
% % infidx = isinf(xgd);
% % xgd(infidx) = [];
% % vds(infidx) = [];
% % 
% % xfun = @(x,xdata)x(1)*sinh(x(2)*(xdata-x(3)))+x(4);
% % % x0 = [.01 .01 250 1];
% % x0 = [10 .01 250 60];
% % 
% % xpar = lsqcurvefit(xfun,x0,vds,xgd);
% % xgd_fit = xfun(xpar,vds);
% % 
% % figure
% % hold on
% % plot(vds, xgd)
% % plot(vds, xgd_fit)
% % hold off
% % 
% % figure
% % plot(vds, cgd(1:end-2));
% % hold on
% % plot(vds, -1/beta_gd*atanh((xgd_fit-gamma_gd)/alpha_gd-1))
% % hold off
