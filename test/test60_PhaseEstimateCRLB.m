%% TEST SCRIPT FOR CALCULATING CRAMER_RAO LOWER BOUND OF PHASE ESTIMATION
% THE CRLB IS THE LOWER BOUND OF MVU ESTIMATION
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
%% Estimate constant phase
% By applying the parameter transformation $\theta = exp(j*\phi)$, the data
% model is linearized. The MVU estimator of transformed parameter can be
% obtained easily. Since it is also an ML estimator, by applying the
% invariance property of ML estimation, the true phase estimation can be
% obtained

clear

bitPerSymbol = 2;
baudRate = 1e9;
sampleRate = 8*baudRate;

RUNN = 500;

% signal power
ps = 2;

snr = -10:10; % in dB

% noise power in dB
sigma2 = dbw(ps)-snr;

theta_hat = ones(RUNN,length(sigma2));

%% Frame Length 1024
frameLength = 1024;

for run=1:RUNN

txBits = randi([0 1],bitPerSymbol,frameLength);
txBaud = symbolizerGrayQam(txBits);

for k=1:length(sigma2)
    
    wns = genWGN(size(txBaud,1),size(txBaud,2),sigma2(k),'dBw','complex');
    
    theta = exp(1j*pi/5);
    
    % signal model
    txx = theta.*txBaud + wns;
    
    % estimator
    theta_hat(run,k) = sum(txx(1:frameLength).*conj(txBaud(1:frameLength))) / sum(txBaud(1:frameLength).*conj(txBaud(1:frameLength)));
end
end

% transform
phi_hat = angle(theta_hat);

% this is CRLB for theta
CRLB = 1./(frameLength.*idbw(snr));
% this is CRLB for phi
CRLBt = 1./(frameLength.*2.*idbw(snr));

figure(1); 
semilogy(snr,CRLB,'-',snr,var(theta_hat),'o'); hold on;
semilogy(snr,CRLBt,'-',snr,var(phi_hat),'d');
xlabel('SNR [dB]'); ylabel('Estimation variance');  grid on

pause(1);

%% Frame Length 2048
frameLength = 4096;

for run=1:RUNN

txBits = randi([0 1],bitPerSymbol,frameLength);
txBaud = symbolizerGrayQam(txBits);

for k=1:length(sigma2)

    wns = genWGN(size(txBaud,1),size(txBaud,2),sigma2(k),'dBw','complex');
    
    theta = exp(1j*pi/5);
    
    % signal model
    txx = theta.*txBaud + wns;
    
    % estimator
    theta_hat(run,k) = sum(txx(1:frameLength).*conj(txBaud(1:frameLength))) / sum(txBaud(1:frameLength).*conj(txBaud(1:frameLength)));
end
end

% transform
phi_hat = angle(theta_hat);

% this is CRLB for theta
CRLB = 1./(frameLength.*idbw(snr));
% this is CRLB for phi
CRLBt = 1./(frameLength.*2.*idbw(snr));

figure(1); 
semilogy(snr,CRLB,'-',snr,var(theta_hat),'s'); hold on;
semilogy(snr,CRLBt,'-',snr,var(phi_hat),'^');
xlabel('SNR [dB]'); ylabel('Estimation variance');  grid on

% larger frame length, lower crlb
