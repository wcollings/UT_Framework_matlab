function [ pos, f, alpha ] = template_pars( sig_t, sig_d )
% TEMPLATE_PARS takes a signal sig_d and a time vector sig_t and computes a
% set of template parameters which describe some of the signal's features. 
%   

% Determing if the isolated waveform is a rising or falling edge. It is
% assumed that the signal provided consists of a single, complete
% transition without anomolous behavior such as self-sustained oscillation.
rise = sig_d(end)-sig_d(1) > 0;

% Locate the index corresponding to the event edge. If no edge is detected,
% use the first value as the event edge. 
sig_d = fillmissing(sig_d, 'nearest');
sw_event = findchangepts(sig_d,'MaxNumChanges',1);
if isempty(sw_event)
    sw_event = 1;
end

% Crop the signal data preceding the switching event. This is necessary to
% isolate the transient portion of the switching event in order to extract
% template parameters in a more targetted fashion. 
sig_t = sig_t(sw_event:end);
sig_d = sig_d(sw_event:end);

% To compute percent overshoot, a steady state offset must be extracted.
% If the signal is known to reach steady state or nearly steady state at
% the end of the isolated signal, the final value can be used. As it is not
% known if this is the case, the mean value of the second half of the
% signal will be used. This allows for oscillation but assumes that the
% ringing frequency is such that the ringing period is much shorter than
% the isolated signal length.
mid = fix(length(sig_d)/2);
ss = mean(sig_d(mid:end));
% ss = sig_d(end);                % Use final value as steady state offset

% If a rising edge is detected, the overshoot amount is the maximum value
% minus the steady state value. If a rising edge is not detected, a falling
% edge is assumed and the steady state value minus the minimum value is
% used. This assumes underdamped, second-order-system-like behavior at the
% switching edge.
if rise
    os = max(sig_d) - ss;
else
    os = ss - min(sig_d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Overshoot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos = os/ss;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the frequency content of the signal using an FFT. Use the
% dominant frequency from the FFT as the ringing frequency. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(sig_d);              % Length of signal
n = 2^nextpow2(L);              % Pad signal
signal = sig_d - ss;            % Center around x axis
%%%%%%%%%%%%%%%%%%%%%%%
[y, idxs] = findpeaks(signal);
x = sig_t(idxs);
if length(x) == 1
    x(2)=sig_t(end);
    y(2)=signal(end);
end
dy = abs(diff(y));
dx = diff(x);
w1i = mean([y(1:end-1);y(2:end)])/mean(y);
w2i = dy/mean(dy);
w3i = (x(end)-x(1))./(mean([x(1:end-1);x(2:end)])-x(1));
w4i = dx/mean(dx);
w = w1i.*w2i.*w3i.*w4i;
f = sum(w)/sum(w.*dx);
%%%%%%%%%%%%%%%%%%%%%%%
% % Y = fft(signal,n);              % Take fft
% % P2 = abs(Y/n);                  % Double sided amplitude spectrum
% % P1 = P2(1:n/2+1);               % Compute single sided amplitude spectrum
% % P1(2:end-1) = 2*P1(2:end-1);    % Double inner values
% % T = sig_t(2) - sig_t(1);        % Sampling period
% % Fs = 1/T;                       % Sampling frequency
% % freqs = 0:(Fs/n):(Fs/2-Fs/n);   % Frequency vector
% % fft_d = P1(1:n/2);              % Single sided spectrum
% % [~,I] = max(fft_d);             % Extract index of dominent frequency
% % %%%%%%%%%%%%%%%%%
% % %%% Frequency %%%
% % f = freqs(I);%%%%
% % if iscolumn(fft_d)
% %     fft_d = fft_d';
% % end
% % f_val = vertcat(fft_d, freqs);
% % %%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%

% Fit the attenuation of the signal. Use the upper peak envelope smoothed
% over the number of samples in two full ringing periods as the curve to
% fit.
sig_T = (sig_t(end) - sig_t(1));
n_T = sig_T*f;
dur_samps = L/n_T; 
samps = fix(2*dur_samps);
[upp, ~] = envelope(sig_d, samps, 'peak');  

% For the exponential curve fitting, assumptions are made about signal
% features. As before, it is assumed that the steady state and overshoot
% values are known. The starting point assumption is that the peak value
% occurs at the first sample. This may be untrue, but is approximately 
% correct if the dominant ringing period is much shorter than the current 
% signal length. It is also assumed, as stated before, that the signal
% approximately reaches steady state behavior. The initial guess for the
% attenuation constant is 1/(0.25*sig_T);


% % % % % alpha_initial = 4/sig_T;
% % % % % kinitial = [ss, os, alpha_initial];
% % % % % opts = optimset('Display', 'off');
% % % % % exp_fit = lsqcurvefit(@(k,m)k(1)+k(2)*exp(-k(3)*m), kinitial, (sig_t-sig_t(1))', upp',[],[], opts);  % Fit exponential curve w/ least squares


%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_initial = 4/sig_T;
kinitial = [ss, os, alpha_initial, x(1)];
opts = optimset('Display', 'off');
% % % exp_fit = lsqcurvefit(@(k,m)k(1)+k(2)*exp(-k(3)*m), kinitial, (sig_t-sig_t(1))', upp',[],[], opts);  % Fit exponential curve w/ least squares
exp_fit = lsqcurvefit(@(k,m)k(1)+k(2)*exp(-k(3)*(m-k(4))), kinitial, x, y,[],[], opts);  % Fit exponential curve w/ least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
%%%%% Alpha %%%%%
alpha = exp_fit(3);%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% disp(exp_fit(1))
% disp(exp_fit(2))


% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % plot([sig_t(1), sig_t(end)], [ss, ss], '-r', 'linewidth', 1.2)
% % % plot([sig_t(1), sig_t(end)], [ss+os, ss+os], '-r', 'linewidth', 1.2)
% % % hold off
% % % 
% % % f_fit = ss+os*sin(2*pi*f*sig_t);
% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % plot(sig_t, f_fit, '-r', 'linewidth', 1.2)
% % % hold off
% % % 
% % % alpha_fit = exp_fit(1) + exp_fit(2)*exp(-exp_fit(3)*(sig_t-sig_t(1)));
% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % % % % plot(sig_t, upp, ':k', 'linewidth', 1.5)
% % % plot(x, y+ss, ':k', 'linewidth', 1.5)
% % % % % % plot(sig_t, alpha_fit, '-r', 'linewidth', 1.2)
% % % plot(sig_t, alpha_fit+ss, '-r', 'linewidth', 1.2)
% % % hold off


% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % plot([sig_t(1), sig_t(end)], [ss, ss], '-r', 'linewidth', 1.2)
% % % plot([sig_t(1), sig_t(end)], [ss+os, ss+os], '-r', 'linewidth', 1.2)
% % % hold off
% % % 
% % % f_fit = ss+os*sin(2*pi*f*sig_t);
% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % plot(sig_t, f_fit, '-r', 'linewidth', 1.2)
% % % hold off
% % % 
% % % alpha_fit = exp_fit(1) + exp_fit(2)*exp(-exp_fit(3)*(sig_t-sig_t(1)));
% % % figure
% % % hold on
% % % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % % grid on
% % % xlim([sig_t(1), sig_t(end)]) 
% % % plot(sig_t, upp, ':k', 'linewidth', 1.5)
% % % plot(sig_t, alpha_fit, '-r', 'linewidth', 1.2)
% % % hold off

% % figure
% % subplot(3, 1, 1)
% % hold on
% % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % xlim([sig_t(1), sig_t(end)]) 
% % plot([sig_t(1), sig_t(end)], [ss, ss], '-r', 'linewidth', 1.2)
% % plot([sig_t(1), sig_t(end)], [ss+os, ss+os], '-r', 'linewidth', 1.2)
% % grid on
% % hold off
% % n_dec = 1; 
% % perc_os = round(os/ss*100, n_dec); 
% % title(['Percent Overshoot = ', num2str(perc_os), '%'])
% % 
% % f_fit = ss+os*sin(2*pi*f*sig_t);
% % subplot(3, 1, 2)
% % hold on
% % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % grid on
% % xlim([sig_t(1), sig_t(end)]) 
% % plot(sig_t, f_fit, '-r', 'linewidth', 1.2)
% % hold off
% % title(['Dominant Frequency = ', num2str(f,'%10.2e\n'), 'Hz'])
% % 
% % alpha_fit = exp_fit(1) + exp_fit(2)*exp(-exp_fit(3)*(sig_t-sig_t(1)));
% % subplot(3, 1, 3)
% % hold on
% % plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% % grid on
% % xlim([sig_t(1), sig_t(end)]) 
% % plot(sig_t, upp, ':k', 'linewidth', 1.5)
% % plot(sig_t, alpha_fit, '-r', 'linewidth', 1.2)
% % hold off
% % title(['\alpha = ', num2str(alpha,'%10.2e\n')])
% 
% figure
% hold on
% plot(sig_t, sig_d, '-k', 'linewidth', 1.2)
% xlim([sig_t(1), sig_t(end)]) 
% plot([sig_t(1), sig_t(end)], [ss+os, ss+os], '-b', 'linewidth', 1.2)
% grid on
% hold off
% n_dec = 1; 
% perc_os = round(os/ss*100, n_dec); 
% % title(['Percent Overshoot = ', num2str(perc_os), '%'])
% 
% f_fit = ss+os*sin(2*pi*f*sig_t);
% hold on
% grid on
% xlim([sig_t(1), sig_t(end)]) 
% plot(sig_t, f_fit, '-g', 'linewidth', 1.2)
% hold off
% % title(['Dominant Frequency = ', num2str(f,'%10.2e\n'), 'Hz'])
% 
% alpha_fit = exp_fit(1) + exp_fit(2)*exp(-exp_fit(3)*(sig_t-sig_t(1)));
% hold on
% % grid on
% xlim([sig_t(1), sig_t(end)]) 
% plot(sig_t, upp, ':k', 'linewidth', 1.5)
% plot(sig_t, alpha_fit, '-r', 'linewidth', 1.2)
% hold off
% % title(['\alpha = ', num2str(alpha,'%10.2e\n')])


end

