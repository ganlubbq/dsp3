% TEST SCRIPT FOR COMPARING VARIOUS CARRIER PHASE ESTIMATION ALGORITHM
%
% - Maximum likelihood
% - Stochastic gradient descent (PLL)
%

% M-QAM WITH TIME VARYING PHASE ERROR
function [ber] = CarrierPhaseRecovery(bitpersym, snr, pnvar, stepsize, blocksize)
if nargin < 1
    bitpersym = 2;
end
if nargin < 2
    snr = 20;
end
if nargin < 3
    pnvar = 1E-5;
end
if nargin < 4
    stepsize = 0.02;
end
if nargin < 5
    blocksize = 20;
end

% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

mn = 2 ^ bitpersym;
symlen = 2 ^ 18;

bitTx = randi([0 1], bitpersym, symlen);
symTx = symbolizerGrayQam(bitTx);
symRef = symTx;

% input power
sp = calcrms(symRef).^2;

% define phase noies
txLaserPnVar = 2 * pi * pnvar;
phaseNoise = genLaserPhaseNoise(symlen, txLaserPnVar, 0);
phaseNoise = phaseNoise(:);
% debug, one could also test fixed phase error
% phaseNoise = 0; 

symTx = symTx .* exp(1i * phaseNoise);

% add wgn
th2 = 10*log10(sp) - snr;
z = wgn(size(symTx,1), size(symTx,2), th2, 'dbw', 'complex');
symTxPn = symTx + z;

theta_ini = 0;
thetaPLL = estimateCarrierPhaseLMS(symRef, symTxPn, stepsize, theta_ini);

thetaML = estimateCarrierPhaseML(symRef, blocksize, mn);

symRecPLL = symTxPn .* exp(-1i * thetaPLL);
symRecML = symTxPn .* exp(-1i * thetaML);

bitRxPLL = slicerGrayQam(normalizeQam(symRecPLL,mn), mn);
bitRxML = slicerGrayQam(normalizeQam(symRecML,mn), mn);

ber(1) = sum(abs(bitTx(:) - bitRxPLL(:))) / (symlen * bitpersym);
ber(2) = sum(abs(bitTx(:) - bitRxML(:))) / (symlen * bitpersym);

varEstErrPLL = calcrms(phaseNoise - thetaPLL(:))^2
varEstErrML = calcrms(phaseNoise - thetaML(:))^2

% bert = T_BER_SNR_mQAM(idbw(snr), mn);
% 
% h1=figure; grid on
% plot(snr, bert, snr, berPLL, snr, berML); 
% xlabel('SNR'); ylabel('BER');
% 
% h2=figure; grid on; hold on
% plot(1:symlen,phaseNoise,'LineWidth',2);
% plot(1:symlen,thetaPLL); 
% plot(1:symlen,thetaML); 
% 
% mngFigureWindow(h1,h2);

return
