function [pos, alpha, f] = template(t,sig)
%%% Receive a second-order like switching waveform and return the
%%% corresponding template parameters. 
rise = sig(end) > sig(1); % Check if signal is rising or falling
cp = findchangepts(sig,'maxnumchanges',1);
m = mean(sig(cp:end));
meanlin = m*ones(size(sig));
if rise
    os = max(sig)-m;
    bpidx = find(sig>=meanlin,1)-2; % Back-projected mean line
else
    os = m-min(sig); 
    bpidx = find(sig<=meanlin,1)-2;
end
pos = 100*abs(os/m);
y = sig(bpidx:end);
x = t(bpidx:end);
ml = m*ones(size(y)); 
s = sign(y-ml);
d = abs(diff(s));
dy = diff(y)./diff(x);
[upp, ~] = envelope(y,round(length(y)/20));   % take the envelope using 5% of the datapoints
expfun = @(x,xdata)x(1)+x(2)*exp(-x(3)*(xdata-x(4)));
x0 = [m, os, sign(rise-.5)*10/(x(end)-x(1)), t(bpidx)];
exppar = lsqcurvefit(expfun, x0, x, upp);
expfit = expfun(exppar, x);
% figure
% hold on
% plot(x,y)
% plot(x,expfit)
% plot(x,upp)
% hold off
didx = find(d); 
diffs = diff(mean([x(didx);x(didx+1)]));    % presumed half-period times
weights = abs(diff(expfit(didx)./max(expfit)));      % associated weights
halfper = sum(diffs.*weights)/sum(weights); % weighted average half period
f = 1/(2*halfper);                          % computed dominant frequency
% figure
% hold on
% plot(x,y)
% plot(x,610+75*sin(x*2*pi*f))
% hold off
% [ya,xa] = findpeaks(smooth(y,11),x);
% sigpow = nlfilter(y,[1,30],@rms);
att = (y-m).^2;
[ya,xa] = findpeaks(att,x);
n = 100; 
sampx = logspace(log10(xa(1)),log10(xa(end)),n);
ya = pchip(xa,ya,sampx);
xa = sampx;
x0 = [0, max(ya), sign(rise-.5)*10/(x(end)-x(1)), xa(1)];
attpar = lsqcurvefit(expfun,x0,xa,ya);
attfit = expfun(attpar,xa);
% figure
% hold on
% plot(xa,ya)
% plot(xa,attfit)
% hold off
alpha = attpar(3); 
end



