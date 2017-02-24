%% TEST SCRIPT FOR CARRIER PHASE ESTIMATION WITH PHASE-LOCK LOOP 
% LOW SNR IS CONSIDERED
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
%% QPSK WITH TIME VARYING PHASE ERROR
clear

% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

bitpersym = 2;
mn = 2^bitpersym;
symlen = 2^16;

% input power
a = constellation(mn);
sp = sum(abs(a).^2) / mn;

%% Low SNR range
snr = -10:0;

% modulation
bitTx = randi([0 1], bitpersym, symlen);
symTx = symbolizerGrayQam(bitTx);
symRef = symTx;

% define phase noies
txLaserPnVar = 2*pi*20e3/2e9;
phaseNoise = genLaserPhaseNoise(symlen, txLaserPnVar, pi/6);
phaseNoise = phaseNoise(:);

% debug, one could also test fixed phase error
% phaseNoise = 0; 

% add phase noise
symTx = symTx .* exp(1i * phaseNoise);

for ii = 1:length(snr)
    
    % add noise
    sigma2 = 10*log10(sp) - snr(ii);
    symTxPn = symTx + genWGN(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');
    
    muV = 10.^(-3.5:0.1:-0.5);
    
    % initialize stochastic gradient descent algorithm, implemented as a
    % phase-lock loop
    for jj = 1:length(muV)
        mu = muV(jj);
        symRec(1) = symTxPn(1);
        phi(1) = 0;
        for k = 2:length(symTxPn)
            % output
            symRec(k) = symTxPn(k) .* exp(-1i * phi(k-1));
            
            % stochastic gradient, also a PED with a typical S-curve
            grad(k) = -imag(symRec(k) .* conj(symRef(k)));
            
            % update filter coeff. along opposite direction of gradient
            phi(k) = phi(k-1) - mu*grad(k);
            
            % squared error
            J(k) = abs(symRec(k) - symRef(k)).^2;
        end
        symRx = normalizeQam(symRec, mn);
        bitrx = slicerGrayQam(symRx, mn);
        ber(jj,ii) = nnz(bitTx-bitrx) / (symlen*2);
        varEstErrL(jj,ii) = calcrms(phaseNoise - phi(:))^2;
    end
end

% bert = T_BER_SNR_mQAM(idbw(snr),mn);
% h1=figure; plot(snr,bert,snr,ber); grid on
% xlabel('SNR'); ylabel('BER');
% h2=figure; plot(1:symlen,phaseNoise,'LineWidth',2); hold on
%            plot(1:symlen,phi); grid on
% mngFigureWindow(h1,h2,'lr');

% plot phase estimation error vs. stepsize
figure; plot(log10(muV),log10(varEstErrL)); grid on
xlabel('Log10 Stepsize'); ylabel('Log10 Variance of phase est. err.');

% phase estimation error floor can be observed and the floor rises as the
% laser noise variance increases

bert = T_BER_SNR_mQAM(idbw(snr), mn);

snrPenalty = calcSnrBerPenalty(snr', bert', ber', mean(bert));
% plot snr penalty vs stepsize
figure; plot(log10(muV),snrPenalty); grid on
xlabel('Log10 Stepsize'); ylabel('SNR Penalty');

% phase estimation error floor transfers into snr penalty, and this is the
% limition of PLL

% EOF
