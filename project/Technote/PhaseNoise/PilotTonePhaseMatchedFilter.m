%% Optimal filter for phase estimation based on pilot-tone
% Phase estimation algorithm based on pilot-tone contains 2 steps: LPF and
% phase rotation. Define the signal model as x = d * exp(i*pn) + w, the
% target here is to get the optimal snr for the data part, i.e., d *
% exp(i*pn), and according to the theory of matched filtering, we should
% use a filter with the same shape as the data part. However, the true
% spectrum shape of data part is unknown and one could adapt the width of
% the LPF towards the maximal dc component of the filtered data part, which
% is NOT PRACTICAL in real cases...

clear
close all

fc = 1e6;
fs = 20e6;
nsample = 10^6;
t = 0 : (1/fs) : (nsample-1)/fs;


%% random walk phase noise
pn = phase_noise(nsample, 2*pi*1E-3, 0);
% a white gaussian additive noise
an = gaussian_noise(1, nsample, .3, 'linear', 'complex');
% data model with zero freq pilot tone
x = exp(1i * pn) + an;


%% test the best width of filter
filtbw = 10 .^ (5 : 0.1 : log10(fs));
for ii = 1 : length(filtbw)
%     % moving average
%     ntaps = 5;
%     taps = ones(1, ntaps) / ntaps;
%     xf = filter(taps, 1, x);

    % raised cosine
    H = frequency_response(nsample, fs, 0.1, filtbw(ii), 'rc');
    xf = ifft(fft(x) .* H.');

    % gaussian snr with signal power of 1
    % one can observe the performance of phase recovery is not related to
    % the gaussian snr.
    anf = ifft(fft(an) .* H.');  
    snr = dbw(1 / calcrms(anf).^2);
    fprintf('SNR is %.2f dB \n', snr);
    
    % remove the phase noise
    xc = x .* conj(sign(xf));
    pc = exp(1i * pn) .* conj(sign(xf));
    
    mx(ii) = abs(mean(xc))^2;
    mp(ii) = abs(mean(pc))^2;
end
figure; plot(log10(filtbw), mx, 'o-'); hold on
plot(log10(filtbw), mp, 's-'); grid on
xlabel('Log frequency'); ylabel('DC power');
legend('data with noise part', 'data part only');


%% adapt using PID loops
filtbw = 1E5;
stepsize = 10;
mp = 0;
err = [];
for ii = 2 : 300
    H = frequency_response(nsample, fs, 0.1, filtbw, 'rc');
    xf = ifft(fft(x) .* H.');
    pc = exp(1i * pn) .* conj(sign(xf));
    % towards the maximal dc power in data part
    mp(ii) = abs(mean(pc))^2;
    err(ii - 1) = abs(mean(pc))^2 - mp(ii - 1);
    filtbw = filtbw + stepsize * err(ii - 1) * 1E5; % + 0.1 * sum(err) * 1E5;
end
figure; plot(err, 'o-'); hold on
plot(mp, 's-'); grid on
    
% pn_est = calcrms(xc - mean(xc)).^2;
% fprintf('estimated power of noise is %.4g \n', pn_est)

% pc = exp(1i*pn) .* conj(xf) ./ abs(xf);
% nc = an .* conj(xf) ./ abs(xf);

% figure(1); spectrumAnalyzer(nc, fs); hold on
% figure(1); spectrumAnalyzer(xc, fs);
% figure(1); spectrumAnalyzer(pc, fs); hold off

% figure(99); clf; hold on
% psd = spectrumAnalyzer(x, fs);
% psd = spectrumAnalyzer(xf, fs); 
% box on; hold off 

% scatterplot(x)
% scatterplot(xc)

% figure(10); clf; hold on
% plot(pn(1:150), 'LineWidth', 2);
% plot(unwrap(angle(xf(1:150)))); 
% legend('Actual phase noise', 'Estimated phase noise');
% grid on; box on;

