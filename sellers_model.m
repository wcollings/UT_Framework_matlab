% %%% Set Up Workspace
% 
% % Add nested directories to path
% folder = fileparts(which(mfilename));
% addpath(genpath(folder));
% 
% % Import Device Data
% device = 'LSIC1M0120E0080';
% load(device)
% 
% lw = 1.2;
% fsz = 10;
% vds = [vds_rev, vds_fwd];
% id_n55 = [id_rev_n55, id_fwd_n55];
% id_25 = [id_rev_25, id_fwd_25];
% id_150 = [id_rev_150, id_fwd_150];
% 
% idlim = [-80,80];
% vdlim = [-7,10];
% 
% figure
% plot(vds, id_n55, 'linewidth', lw)
% xlim(vdlim)
% ylim(idlim)
% title('Conduction T = -55{\circ}C')
% xlabel('V_D_S (V)')
% ylabel('I_D (A)')
% grid on
% box on
% set(gca, 'fontweight', 'bold', 'fontsize', fsz)
% 
% figure
% plot(vds, id_25, 'linewidth', lw)
% xlim(vdlim)
% ylim(idlim)
% title('Conduction T = 25{\circ}C')
% xlabel('V_D_S (V)')
% ylabel('I_D (A)')
% grid on
% box on
% set(gca, 'fontweight', 'bold', 'fontsize', fsz)
% 
% figure
% plot(vds, id_150, 'linewidth', lw)
% xlim(vdlim)
% ylim(idlim)
% title('Conduction T = 150{\circ}C')
% xlabel('V_D_S (V)')
% ylabel('I_D (A)')
% grid on
% box on
% set(gca, 'fontweight', 'bold', 'fontsize', fsz)
% 
% % % figure
% % % plot(vds, r_n55, 'linewidth', lw)
% % % xlim(vdlim)
% % % ylim(idlim)
% % % title('Resistance T = -55{\circ}C')
% % % xlabel('V_D_S (V)')
% % % ylabel('R_D (\Omega)')
% % % grid on
% % % box on
% % % set(gca, 'fontweight', 'bold', 'fontsize', fsz)
% 
% % % id = flipud(id_fwd_n55);
% % id = flipud(id_fwd_25);
% % % id = flipud(id_fwd_150);
% % close all
% % model_coefficients = sellers_forward(vds_fwd, vgs_fwd, id);
% 
% 
% % r_n55 = vds./id_n55;
% % r_25 = vds./id_25;
% % r_150 = vds./id_150;
% % figure
% % plot(vds, r_n55, 'linewidth', lw)
% % figure
% % plot(vds, r_25, 'linewidth', lw)
% % figure
% % plot(vds, r_150, 'linewidth', lw)
%  
% % dv = mean([vds(1:end-1);vds(2:end)]);
% % diff_v = diff(vds);
% % diff_id_n55 = diff(id_n55,1,2);
% % diff_id_25 = diff(id_25,1,2);
% % diff_id_150 = diff(id_150,1,2);
% % r_n55 = diff_v./diff_id_n55;
% % r_25 = diff_v./diff_id_25;
% % r_150 = diff_v./diff_id_150;
% % figure
% % plot(dv, r_n55, 'linewidth', lw)
% % figure
% % plot(dv, r_25, 'linewidth', lw)
% % figure
% % plot(dv, r_150, 'linewidth', lw)
% 
% 
% % %%%%%%%%%%%%%%%%
% % close all
% % clear
% % x = 0:0.1:100; 
% % R = 1; 
% % y1 = x./R; 
% % figure
% % plot(x,y1,'-k','linewidth',1.2)
% % 
% % y2 = y1 + 20; 
% % 
% % hold on
% % plot(x,y2,'--k','linewidth',1.2)
% % hold off
% % xlabel('V')
% % ylabel('I')
% % grid on
% % box on
% % set(gca, 'fontweight', 'bold')
% % 
% % R1 = y1./x;
% % R2 = y2./x;
% % figure
% % hold on
% % plot(x, R1,'-k','linewidth',1.2)
% % plot(x, R2,'--k','linewidth',1.2)
% % hold off
% % xlabel('V')
% % ylabel('R')
% % grid on
% % box on
% % set(gca, 'fontweight', 'bold')
% 
% 
% %%%%%%%%%%% Reverse Conduction
% % Rfun = @(x,xdata)x(1)*(1/x(2)*log(cosh(x(2)*(x(3)-xdata)))+xdata)+x(4);
% % % Rfun = @(x,xdata)x(1)/(2*x(2))*(-log(x(2)^2*(x(3)-xdata).^2+1)+2*x(2).*(x(3)-xdata).*atan(x(2)*(x(3)-xdata))+pi*x(2)*xdata)+x(4);
% % % Rfun = @(x,xdata)x(1)*(log(cosh(x(2)*(x(3)-xdata)))+x(2)*xdata)/x(2);
% % % Sfun = @(x,xdata)x(1)*(tanh(x(2)*(xdata-x(3)))+1);
% % 
% % x0 = [1, 1, 1, 0];
% % % x0 = [10    .5    5  -50];
% % 
% % x1 = fliplr(-vds_rev);
% % y1 = fliplr(-id_rev_25(4,:));
% % % y1 = fliplr(-id_rev_25(6,:));
% % 
% % Rpar1 = lsqcurvefit(Rfun,x0,x1,y1);
% % 
% % % x = -10:0.01:10;
% % % Rsta1 = Rfun(x0,x);
% % % Ssta1 = Sfun(x0,x);
% % % Rfit1 = Rfun(Rpar1,x);
% % 
% % Rsta1 = Rfun(x0,x1);
% % % Ssta1 = Sfun(x0,x1);
% % Rfit1 = Rfun(Rpar1,x1);
% % 
% % figure
% % hold on
% % plot(x1,y1)
% % % plot(x1,Rsta1)
% % % plot(x, cumtrapz(x,Ssta1))
% % plot(x1,Rfit1)
% % % plot(x,Ssta1)
% % hold off
% % 
% % % figure
% % % plot(x,Ssta1)
% % 
% % 
% % dx = mean([x1(1:end-1);x1(2:end)]);
% % dy = diff(y1)./diff(x1);
% % 
% % figure
% % plot(dx,dy)
% 
% % Fdiode = @(x,xdata)x(1)*(exp(x(2)*xdata)-1);
% % 
% % x0 = [1, .026];
% % 
% % x = fliplr(-vds_rev);
% % y1 = fliplr(-id_rev_25(6,:));
% % 
% % Dpar = lsqcurvefit(Fdiode, x0, x, y1);
% % Dfit = Fdiode(Dpar,x);
% % 
% % figure
% % hold on
% % plot(x,y1)
% % plot(x,Dfit)
% % hold off
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LTSPICE compact diode model: 
% A diode requires a .model card to specify its characteristics. There are 
% two types of diodes available. One is a conduction region-wise linear 
% model that yields a computationally light weight representation of an 
% idealized diode. It has three linear regions of conduction: on, off and 
% reverse breakdown. Forward conduction and reverse breakdown can 
% non-linear by specifying a current limit with Ilimit(revIlimit). tanh() 
% is used to fit the slope of the forward conduction to the limit current. 
% The parameters epsilon and revepsilon can be specified to smoothly switch
% between the off and conducting states. A quadratic function is fit 
% between the off and on state such that the diode's IV curve is continuous
% in value and slope and the transition occurs over a voltage specified by
% the value of epsilon for the off to forward conduction and revepsilon for
% the transition between off and reverse breakdown.

%%% Model parameters: 
% Ron = 1
% Roff = 1/Gmin
% Vfwd = 0
% Vrev = Infin
% Rrev = Ron
% Ilimit = Infin
% Revilimit = Infin
% Epsilon = 0
% Revepsilon = 0

% epsilon = 1;
% Ron = 1;
% 
% a = 1/(2*epsilon*Ron);
% % fdiode = @(xdata)(xdata>=0&&xdata<epsilon).*a*xdata.^2+(xdata>epsilon).*1/Ron*xdata;
% 
% % fdiode = @(xdata)piecewise(xdata<=0,0,xdata>0&xdata<=epsilon,a*xdata.^2,xdata>epsilon,xdata/Ron);
% fdiode = @(v)((v>0)&(v<epsilon)).*(a*v.^2)+(v>=epsilon).*(v/Ron-a*epsilon^2);
% 
% syms('t');
% y(t) = piecewise(t<=0,0,t>0&t<=epsilon,a*t.^2,t>epsilon,t/Ron-a*epsilon^2);
% 
% 
% Add nested directories to path
folder = fileparts(which(mfilename));
addpath(genpath(folder));

% Import Device Data
device = 'LSIC1M0120E0080';
load(device)

x = fliplr(-vds_rev);
% y1 = fliplr(-id_rev_25(6,:));

% 
% rawdata = LTspice2Matlab('DiodeTest');
% headings = rawdata.variable_name_list; 
% idiode_idx = strcmpi(headings, 'Ix(d:A)');
% idiode = rawdata.variable_mat(idiode_idx,:);
% v = rawdata.time_vect;
% 
% figure
% hold on
% plot(v, idiode)
% plot(v, y(v))
% plot(v, fdiode(v))
% hold off

% a = (1/(2*x(1)*x(2)));
% fdiode = @(x,v)((v<=0).*0+(v>0)&(v<x(1))).*((1/(2*x(1)*x(2)))*v.^2)+(v>=x(1)).*(v/x(2)-(1/(2*x(1)*x(2)))*x(1)^2);
fdiode = @(x,v)(((v-x(3))<=0).*0+((v-x(3))>0)&((v-x(3))<x(1))).*((1/(2*x(1)*x(2)))*(v-x(3)).^2)+((v-x(3))>=x(1)).*((v-x(3))/x(2)-(1/(2*x(1)*x(2)))*x(1)^2);
epsilon0 = 5;
Ron0 = 50e-3;
Vfwd0 = 3;
x0 = [epsilon0, Ron0, Vfwd0];
lb = [1e-6 1e-6 -inf]; 
ub = [inf inf inf];
dpar = lsqcurvefit(fdiode,x0,x,y1,lb,ub);
dfit = fdiode(dpar,x);

% rawdata = LTspice2Matlab('DiodeTest');
% headings = rawdata.variable_name_list; 
% idiode_idx = strcmpi(headings, 'Ix(d:A)');
% idiode = rawdata.variable_mat(idiode_idx,:);
% v = rawdata.time_vect;

figure
hold on
plot(x,y1)
plot(x,dfit)
% plot(v,idiode)
% plot(x,fdiode(x0,x))
hold off

