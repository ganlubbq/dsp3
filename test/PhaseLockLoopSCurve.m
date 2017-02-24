%% TEST SCRIPT FOR PHASE-LOCKED LOOP S-CURVE
% Observe the PLL S-Curve under different SNR. The instanstanuous gradient
% gets noisier for lower SNR and requires smaller stepsize (smaller loop
% filter bandwidth) to lower the steady state error
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

%% Observe averaged Phase-locked loop S-curve for various SNR
snr = -20:20;

% define phase error
phaseNoise = -pi:0.1:pi+0.1;

for ii = 1:length(snr)
    
    % generate wgn
    sigma2 = 10*log10(sp) - snr(ii);
    z = wgn(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');
    
    for k = 1:length(phaseNoise)
        % data model
        symTxPn = symTx .* exp(1i * phaseNoise(k)) + z;
        
        % average over random data; normalized wrt the signal power
        sc(k,ii) = (1/sp) * mean(imag(symTxPn .* conj(symRef)));
    end
end

h1=figure; plot(phaseNoise/pi,sc); grid on; xlim([-1 1]);
xlabel('Phase Error (\phi-\phi_k)/\pi'); ylabel('PED mean');

sc = [];

%% Observe instantanuous Phase-locked loop S-curve for low SNR
snr = -10;

for ii = 1:10
    
    % randomly pick an symbol index
    tmp = round(rand()*symlen);
    
    sigma2 = 10*log10(sp) - snr;
    z = wgn(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');
    
    for k = 1:length(phaseNoise)
        
        % data model
        symTxPn = symTx .* exp(1i * phaseNoise(k)) + z;
        
        % use latest data only; normalized wrt the signal power
        sc(k,ii) = (1/sp) * (imag(symTxPn(tmp) .* conj(symRef(tmp))));
    end
end

h2=figure; plot(phaseNoise,sc); grid on; xlim([-pi pi]);
xlabel('Phase Error'); ylabel('PED');

mngFigureWindow(h1,h2);

