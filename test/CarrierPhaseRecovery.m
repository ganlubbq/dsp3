%% TEST SCRIPT FOR CARRIER PHASE ESTIMATION COMPARISON
% not finished

% M-QAM WITH TIME VARYING PHASE ERROR
function [] = CarrierPhaseRecovery(bitpersym)
if nargin < 1
    bitpersym = 2;
end

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

mn = 2 ^ bitpersym;
symlen = 2 ^ 16;

bitTx = randi([0 1], bitpersym, symlen);
symTx = symbolizerGrayQam(bitTx);
symRef = symTx;

% input power
sp = calcrms(symRef).^2;

% define phase noies
txLaserPnVar = 2 * pi * 10e3 / 2e9;
phaseNoise = genLaserPhaseNoise(symlen, txLaserPnVar, 0);
phaseNoise = phaseNoise(:);
% debug, one could also test fixed phase error
% phaseNoise = 0; 

symTx = symTx .* exp(1i * phaseNoise);

snr = -20 : 0;
for ii = 1:length(snr)
    % add wgn
    th2 = 10*log10(sp) - snr(ii);
    z = wgn(size(symTx,1), size(symTx,2), th2, 'dbw', 'complex');
    symTxPn = symTx + z;
    
    stepsize = .05;
    theta_ini = 0;
    thetaPLL = estimateCarrierPhaseLMS(symRef, symTxPn, stepsize, theta_ini);
    
    blocksize = 2;
    thetaML = estimateCarrierPhaseML(symRef, symTxPn, blocksize);
    
    symRecPLL = symTxPn .* exp(-1i * thetaPLL);
    symRecML = symTxPn .* exp(-1i * thetaML);
    
    bitRxPLL = slicerGrayQam(normalizeQam(symRecPLL,mn), mn);
    bitRxML = slicerGrayQam(normalizeQam(symRecML,mn), mn);
    
    berPLL(ii) = nnz(bitTx - bitRxPLL) / (symlen * bitpersym);
    berML(ii) = nnz(bitTx - bitRxML) / (symlen * bitpersym);
end

varEstErrPLL = calcrms(phaseNoise - thetaPLL(:))^2
varEstErrML = calcrms(phaseNoise - thetaML(:))^2

bert = T_BER_SNR_mQAM(idbw(snr), mn);

h1=figure; grid on
plot(snr, bert, snr, berPLL, snr, berML); 
xlabel('SNR'); ylabel('BER');

h2=figure; grid on; hold on
plot(1:symlen,phaseNoise,'LineWidth',2);
plot(1:symlen,thetaPLL); 
plot(1:symlen,thetaML); 

mngFigureWindow(h1,h2);

return
