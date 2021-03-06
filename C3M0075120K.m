%%% Breakdown of the model for the Wolfspeed SiC power MOSFET C3M0075120K



% Since there is temperature dependence of multiple values, set a constant
% T = 25C for analysis
Tj = 25; 



%%% Conduction Behavior:

% % Forward Conduction
% G1 d1 s value {
% +	if(V(s,d3)<0,
% +		0
% +		,
% +		if (V(gk,s)<V(NET3) ,
% +		-((0.035)*(v(gk,s)-V(NET3)))*(-(1+p10*v(s,d3))*0.008)*(((log(1+exp(v(gk,s)-V(NET3))))**2)-
% +		((log(1+exp(v(gk,s)-V(NET3)-(0.854*v(s,d3)))))**2))
% +		,
% +		-((v(NET5)+v(NET4))*(v(gk,s)-V(NET3)))*(1+v(P8)*v(s,d3))*(((log(1+exp(v(gk,s)-V(NET3))))**2)-
% +		((log(1+exp(v(gk,s)-V(NET3)-(V(NET2)*v(s,d3)*((1+exp(-v(NET10)*v(s,d3)))**v(NET1))))))**2))
% +		)
% +			)
% +			}

vds = 0:0.1:20;

vgs = (0:1:20)';
% vgs = (0:.1:20)';


p10  = 0.0011;
p11 = -8;
p12 = 19;
p13 = 15;

vgk = min(max(p11,vgs),p12);

v3 = 8e-6*Tj^2-4.7e-3*Tj+2.8224;

% e5		NET5   0	Value {
% +					if (V(gk,s)>11 ,
% *+				(-0.9267*V(gk,s)**3+49.313*V(gk,s)**2-877.727*V(gk,s)+5351.268)/10000
% +				((87.641n*V(Tj)**3+46.001u*V(Tj)**2-15.03m*V(Tj)-0.13539)*V(gk,s)**3+
% +				(-4.7725u*V(Tj)**3-2.0118m*V(Tj)**2+0.69225*V(Tj)+13.826)*V(gk,s)**2+
% +				(86.016u*V(Tj)**3+27.876m*V(Tj)**2-10.322*V(Tj)-366.33)*v(gk,s)+
% +				(-496.9u*V(Tj)**3-0.12272*V(Tj)**2+49.167*V(Tj)+3084.3))/10000
% +					,
% +					if (V(gk,s)<=11 & V(gk,s)>9,
% *+				(15*V(gk,s)**2-245*V(gk,s)+1470)/10000
% +				((8.3091u*V(Tj)**3+1.2517m*V(Tj)**2-0.30635*V(Tj)-4.25318)*(V(gk,s)**2)+
% +				(-166.98u*V(Tj)**3-21.874m*V(Tj)**2+4.7236*V(Tj)+53.187)*v(gk,s)+
% +				(821.34u*V(Tj)**3+90.986m*V(Tj)**2-15.564*V(Tj)+475.4))/10000
% +					,
% *+				(13.958*V(gk,s)**2-158.333*V(gk,s)+774.375)/10000
% +				((-8.337u*V(Tj)**3+1.507m*V(Tj)**2-94.69m*V(Tj)+2.806)*(V(gk,s)**2)+
% +				(92.64u*V(Tj)**3-16.24m*V(Tj)**2+0.5932*V(Tj)-61.12)*v(gk,s)+
% +				(-166.9u*V(Tj)**3+19.64m*V(Tj)**2+4.464*V(Tj)+932.7))/10000
% +			)
% +			)
% +			}

% % If vgk > 11:
% v5 = ((87.641e-9*Tj^3+46.001e-6*Tj^2-15.03e-3*Tj-0.13539)*vgk.^3+(-4.7725e-6*Tj^3-2.0118e-3*Tj^2+0.69225*Tj+13.826)*vgk.^2+...
%     (86.016e-6*Tj^3+27.876e-3*Tj^2-10.322*Tj-366.33)*vgk+(-496.9e-6*Tj^3-0.12272*Tj^2+49.167*Tj+3084.3))/10000;
% v5 = (-0.4810*vgk.^3+29.8003*vgk.^2+(-605.6135)*vgk+(4.2290e+03))/10000;
% v5 = @(v)(-0.4810*v.^3+29.8003*v.^2+(-605.6135)*v+(4.2290e+03))/10000;

% % If vgk <= 11 & vgk > 9
% v5 = ((8.3091e-6*Tj^3+1.2517e-3*Tj^2-0.30635*Tj-4.25318)*(vgk.^2)+(-166.98e-6*Tj.^3-21.874e-3*Tj^2+4.7236*Tj+53.187)*vgk+...
%     (821.34e-6*Tj^3+90.986e-3*Tj^2-15.564*Tj+475.4))/10000;
% v5 = (-10.9998*vgk.^2+154.9967*vgk+155.9997)/10000;

% % If vgk <= 9
% v5 = ((-8.337e-6*Tj^3+1.507e-3*Tj^2-94.69e-3*Tj+2.806)*(vgk.^2)+(92.64e-6*Tj^3-16.24e-3*Tj^2+0.5932*Tj-61.12)*vgk+...
%     (-166.9e-6*Tj^3+19.64e-3*Tj^2+4.464*Tj+932.7))/10000;
% v5 = (1.2504*vgk.^2+(-54.9925)*vgk+1.0540e+03)/10000;

s5_1 = vgk > 11;
s5_2 = vgk <= 11 & vgk > 9;
s5_3 = vgs <= 9;

v5_1 = s5_1.*((-0.4810*vgk.^3+29.8003*vgk.^2+(-605.6135)*vgk+(4.2290e+03))/10000);
v5_2 = s5_2.*((-10.9998*vgk.^2+154.9967*vgk+155.9997)/10000);
v5_3 = s5_3.*((1.2504*vgk.^2+(-54.9925)*vgk+1.0540e+03)/10000);

v5 = v5_1 + v5_2 + v5_3;

% figure
% plot(vgs, v5)



% e_p8	P8	0	Value {Limit(((95.93n*V(Tj)**3-17.89u*V(Tj)**2+8.478u*V(Tj)+35.59m)*(V(gk,s)**3)+
% +				(-4.135u*V(Tj)**3+831u*V(Tj)**2-8.584m*V(Tj)-2.647)*(V(gk,s)**2)+
% +				(54.29u*V(Tj)**3-11.48m*V(Tj)**2+0.1753*V(Tj)+51.4)*v(gk,s)+
% +				(-216u*V(Tj)**3+46.84m*V(Tj)**2-0.7812*V(Tj)-210.8))/1000,0.001,0.2)
% +					}

% vp8_sig = ((95.93e-9*Tj^3-17.89e-6*Tj^2+8.478e-6*Tj+35.59e-3)*(vgk.^3)+(-4.135e-6*Tj^3+831e-6*Tj^2-8.584e-3*Tj-2.647)*(vgk.^2)+...
%     (54.29e-6*Tj^3-11.48e-3*Tj^2+0.1753*Tj+51.4)*vgk+(-216e-6*Tj^3+46.84e-3*Tj^2-0.7812*Tj-210.8))/1000;
vp8_sig = (0.0261*(vgk.^3)+(-4.135e-6*Tj^3+831e-6*Tj^2-8.584e-3*Tj-2.647)*(vgk.^2)+49.4558*vgk+(-204.4300))/1000;
vp8 = min(max(vp8_sig,0.001),0.2);

v2 = 15.35e-3*vgk+371.85e-3;

% e10		NET10   0	Value {Limit(((-94.87u*V(Tj)**2+25.49m*V(Tj)-0.8726)*(V(gk,s)**3)+
% +					(3.038m*V(Tj)**2-0.8788*V(Tj)+35.82)*(V(gk,s)**2)+
% +					(-29.94m*V(Tj)**2+9.729*V(Tj)-501.7)*v(gk,s)+
% +					(85.19m*V(Tj)**2-34.18*V(Tj)+2452))/1000,0.001,3.7)
% +				      }
% v10_sig = ((-94.87e-6*Tj^2+25.49e-3*Tj-0.8726)*(vgk.^3)+(3.038e-3*Tj^2-0.8788*Tj+35.82)*(vgk.^2)+...
%     (-29.94e-3*Tj^2+9.729*Tj-501.7)*vgk+(85.19e-3*Tj^2-34.18*Tj+2452))/1000;
v10_sig = ((-0.2946)*vgk.^3+15.7488*vgk.^2+(-277.1875)*vgk+1.6507e+03)/1000;
v10 = min(max(v10_sig,0.001),3.7);



% e1		NET1	0	Value {Limit(((997.8n*V(Tj)**3-167.2u*V(Tj)**2+2.679m*V(Tj)-97.69m)*V(gk,s)**4+
% +				(-48.54u*V(Tj)**3+7.754m*V(Tj)**2-69.85m*V(Tj)+2.697)*V(gk,s)**3+
% +				(839.1u*V(Tj)**3-0.1254*V(Tj)**2-0.1785*V(Tj)+1.003)*V(gk,s)**2+
% +				(-6.049m*V(Tj)**3+0.8158*V(Tj)**2+15.49*V(Tj)-400.7)*v(gk,s)+
% +				(15.16m*V(Tj)**3-1.738*V(Tj)**2-88.88*V(Tj)+3393))/1000,0.01,15)
% +					}

% v1_sig = ((997.8e-9*Tj^3-167.2e-6*Tj^2+2.679e-3*Tj-97.69e-3)*vgk.^4+(-48.54e-6*Tj^3+7.754e-3*Tj^2-69.85e-3*Tj+2.697)*vgk.^3+...
%     (839.1e-6*Tj^3-0.1254*Tj^2-0.1785*Tj+1.003)*vgk.^2+(-6.049e-3*Tj^3+0.8158*Tj^2+15.49*Tj-400.7)*vgk+...
%     (15.16e-3*Tj^3-1.738*Tj^2-88.88*Tj+3393))/1000;
v1_sig = ((-0.1196)*vgk.^4+5.0386*vgk.^3+(-68.7236)*vgk.^2+401.9094*vgk+321.6250)/1000;
v1 = min(max(v1_sig,0.01),15);

% Id forward when Vgk < v3
Id_forward_1 = ((0.035)*vgk-v3)*(-(1+p10*vds)*0.008).*(((log(1+exp(vgk-v3))).^2)-((log(1+exp(vgk-v3-(0.854*vds)))).^2));
% Id_forward_1(vgk>=v3,:)=NaN;

% Id forward when Vgk >= v3
Id_forward_2 = ((v5.*(vgk-v3))).*(1+vp8.*vds).*(((log(1+exp(vgk-v3))).^2)-((log(1+exp(vgk-v3-(v2.*vds.*((1+exp(-v10.*vds)).^v1))))).^2));

Id_forward = Id_forward_1.*(vgk<v3)+Id_forward_2.*(vgk>=v3);
% Id_forward_1(vgk>=v3,:)=NaN;

figure
plot(vds, Id_forward)


% % Reverse Conduction
% G2 d1 s value {
% +	if(V(d3,s)<0,
% +		0
% +		,
% +		if (V(gk,s)<V(NET3) ,
% +		((0.035)*(v(gk,s)-V(NET3)))*(-(1+p10*v(d3,s))*0.008)*(((log(1+exp(v(gk,s)-V(NET3))))**2)-
% +		((log(1+exp(v(gk,s)-V(NET3)-(0.854*v(d3,s)))))**2))
% +		,
% +		((v(NET5)*(v(gk,s)-V(NET3))))*(1+v(P8)*v(d3,s))*(((log(1+exp(v(gk,s)-V(NET3))))**2)-
% +		((log(1+exp(v(gk,s)-V(NET3)-(V(NET2)*v(d3,s)*((1+exp(-v(NET10)*v(d3,s)))**v(NET1))))))**2))
% +		)
% +			)
% +			}
