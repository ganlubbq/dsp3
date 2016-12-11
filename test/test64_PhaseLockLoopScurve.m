%% TEST SCRIPT FOR PHASE-LOCKED LOOP S-CURVE
%
%
%% QPSK WITH TIME VARYING PHASE ERROR
clear

% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

bitpersym = 2;
mn = 2^bitpersym;
symlen = 2^14;

% input power
sp = 2;

% modulation
bitTx = randi([0 1],bitpersym,symlen);
symTx = symbolizerGrayQam(bitTx);
symRef = symTx;


%% High SNR Phase-locked loop S-curve
snr = 0:10;

% define phase error
phaseNoise = -pi:0.1:pi+0.1;

for ii = 1:length(snr)
    
    % generate wgn
    sigma2 = 10*log10(sp) - snr(ii);
    z = wgn(1,symlen,sigma2,'dbw','complex');
    
    for k = 1:length(phaseNoise)
        
        % data model
        symTxPn = symTx.*exp(1j*phaseNoise(k))+z;
        
        % average over random data; normalized wrt the signal power
        sc(k,ii) = (1/sp)*mean(imag(symTxPn.*conj(symRef)));
    end
end

figure; plot(phaseNoise/pi,sc); grid on; xlim([-1 1]);
xlabel('Phase Error (\phi-\phi_k)/\pi'); ylabel('PED mean');

sc = [];


%% Low SNR Phase-locked loop S-curve
snr = -10:0;

for ii = 1:length(snr)
    sigma2 = 10*log10(sp) - snr(ii);
    z = wgn(1,symlen,sigma2,'dbw','complex');
    for k = 1:length(phaseNoise)
        % data model
        symTxPn = symTx.*exp(1j*phaseNoise(k))+z;
        % average over random data; normalized wrt the signal power
        sc(k,ii) = (1/sp)*mean(imag(symTxPn.*conj(symRef)));
    end
end

figure; plot(phaseNoise,sc); grid on; xlim([-pi pi]);
xlabel('Phase Error'); ylabel('PED mean');

% the averaged ped output forms sinusoidal s-curve from -pi to pi, and
% holds for both high and low SNR

sc = [];


%% Low SNR Phase-locked loop S-curve instantanuous
snr = -2;

for ii = 1:10
    
    % randomly pick an symbol index
    tmp = round(rand()*symlen);
    
    sigma2 = 10*log10(sp) - snr;
    z = wgn(1,symlen,sigma2,'dbw','complex');
    
    for k = 1:length(phaseNoise)
        
        % data model
        symTxPn = symTx.*exp(1j*phaseNoise(k))+z;
        
        % use latest data only; normalized wrt the signal power
        sc(k,ii) = (1/sp)*(imag(symTxPn(tmp).*conj(symRef(tmp))));
    end
end

figure; plot(phaseNoise,sc); grid on; xlim([-pi pi]);
xlabel('Phase Error'); ylabel('PED');

% the instanstanuous gradient gets noisier for lower snr and requires
% smaller stepsize (smaller loop filter bandwidth) to lower the steady
% state error


% EOF
