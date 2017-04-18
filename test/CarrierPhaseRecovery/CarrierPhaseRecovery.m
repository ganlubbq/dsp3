% TEST SCRIPT FOR COMPARING VARIOUS CARRIER PHASE ESTIMATION ALGORITHM
%
% - Maximum likelihood (ideal for m-psk)
% - Stochastic gradient descent (PLL)
% - Adadelta (variable stepsize PLL)
%
% Adadelta may work well on some irregular cost surface, however not well
% on the time-varying cost surface like the carrier phase noise. Will test
% methods such as Kalman filter.
%
% Test cases:
%   CarrierPhaseRecovery(2, 7, 1E-4, 0.40, 1, 512, 8, 0)
%   CarrierPhaseRecovery(2, 7, 1E-4, 0.04, 1, 512, 8, 1)
%   CarrierPhaseRecovery(4, 10, 1E-4, 0.05, 1, 512, 8, 0)
%   CarrierPhaseRecovery(4, 10, 1E-4, 0.01, 1, 512, 8, 1)
%
% M-QAM WITH TIME VARYING PHASE ERROR
function [ber] = CarrierPhaseRecovery(bitpersym, snr, pnvar, stepsize, blocksize, framesize, trainingsize, ddmode)
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
if nargin < 6
    framesize = 512;
end
if nargin < 7
    trainingsize = 8;
end
if nargin < 8
    ddmode = 0;
end

% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing signals with phase noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pure training mode, all the signals are known
% the maximum likelihood (with block size of 1) is the best
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ddmode == 0
    theta_ini = 0;
    thetaPLL = estimateCarrierPhaseDDPLL(symRef, symTxPn, mn, stepsize, numel(symTxPn), numel(symTxPn), theta_ini);
    thetaAda = estimateCarrierPhaseAdadelta(symRef, symTxPn, mn, numel(symTxPn), numel(symTxPn), theta_ini);
    thetaML = estimateCarrierPhaseML(symRef, symTxPn, blocksize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% training + decision mode, periodic training symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ddmode == 1
    theta_ini = 0;
    thetaPLL = estimateCarrierPhaseDDPLL(symRef, symTxPn, mn, stepsize, framesize, trainingsize, theta_ini);
    thetaAda = estimateCarrierPhaseAdadelta(symRef, symTxPn, mn, framesize, trainingsize, theta_ini);
    thetaML = estimateCarrierPhaseML(symRef, symTxPn, blocksize);
else 
    keyboard;
end

symRecPLL = symTxPn .* exp(-1i * thetaPLL);
symRecAda = symTxPn .* exp(-1i * thetaAda);
symRecML = symTxPn .* exp(-1i * thetaML);

bitRxPLL = slicerGrayQam(normalizeQam(symRecPLL,mn), mn);
bitRxAda = slicerGrayQam(normalizeQam(symRecAda,mn), mn);
bitRxML = slicerGrayQam(normalizeQam(symRecML,mn), mn);

ber(1) = sum(abs(bitTx(:) - bitRxPLL(:))) / (symlen * bitpersym);
ber(2) = sum(abs(bitTx(:) - bitRxAda(:))) / (symlen * bitpersym);
ber(3) = sum(abs(bitTx(:) - bitRxML(:))) / (symlen * bitpersym);

return
