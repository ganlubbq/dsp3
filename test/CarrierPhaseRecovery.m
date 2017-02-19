%% TEST SCRIPT FOR CARRIER PHASE ESTIMATION COMPARISON
% not finished
%
%% QPSK WITH TIME VARYING PHASE ERROR
clear
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
% qpsk
bitpersym = 2;
mn = 2^bitpersym;
% symbol length
symlen = 2^16;
% input power
sp = 2;


%% High SNR Phase-locked loop
snr = 0:10;

bitTx = randi([0 1],bitpersym,symlen);
symTx = symbolizerGrayQam(bitTx);
symRef = symTx;
% define phase noies
txLaserPnVar = 2*pi*10e3/2e9;
phaseNoise = genLaserPhaseNoise(symlen,txLaserPnVar,pi/6);
% phaseNoise = 0; % debug, one could also test fixed phase error

symTx = symTx .* exp( 1j * phaseNoise);

for ii = 1:length(snr)
    th2 = 10*log10(sp) - snr(ii);
    z = wgn(1,symlen,th2,'dbw','complex');
    symTxPn = symTx + z;
    % initialize stochastic gradient descent algorithm, implemented as a
    % phase-lock loop
    mu = 0.01;
    symRec(1) = symTxPn(1);
    phi(1) = 0;
    for k = 2:length(symTxPn)
        % output
        symRec(k) = symTxPn(k).*exp(-1j*phi(k-1));
        % stochastic gradient, also a PED with a typical S-curve
        grad(k) = -imag(symRec(k).*conj(symRef(k)));
        % update filter coeff. along opposite direction of gradient
        phi(k) = phi(k-1) - mu*grad(k);
        % squared error
        J(k) = abs(symRec(k)-symTx(k)).^2;
    end
    symrx = normalizeQam(symRec,mn);
    bitrx = slicerGrayQam(symrx,mn);
    ber(ii) = nnz(bitTx-bitrx)/(symlen*2);
end

varEstErrH = calcrms(phaseNoise-phi)^2

bert = T_BER_SNR_mQAM(idbw(snr),mn);
h1=figure; plot(snr,bert,snr,ber); grid on
xlabel('SNR'); ylabel('BER');
h2=figure; plot(1:symlen,phaseNoise,'LineWidth',2); hold on
           plot(1:symlen,phi); grid on

mngFigureWindow(h1,h2);

%% Low SNR Phase-locked loop
snr = -10:0;

for ii = 1:length(snr)
    th2 = 10*log10(sp) - snr(ii);
    z = wgn(1,symlen,th2,'dbw','complex');
    symTxPn = symTx + z;
    % initialize stochastic gradient descent algorithm, implemented as a
    % phase-lock loop
    mu = 0.001;
    symRec(1) = symTxPn(1);
    phi(1) = 0;
    for k = 2:length(symTxPn)
        % output
        symRec(k) = symTxPn(k).*exp(-1j*phi(k-1));
        % stochastic gradient, also a PED with a typical S-curve
        grad(k) = -imag(symRec(k).*conj(symRef(k)));
        % update filter coeff. along opposite direction of gradient
        phi(k) = phi(k-1) - mu*grad(k);
        % squared error
        J(k) = abs(symRec(k)-symTx(k)).^2;
    end
    symrx = normalizeQam(symRec,mn);
    bitrx = slicerGrayQam(symrx,mn);
    ber(ii) = nnz(bitTx-bitrx)/(symlen*2);
end

varEstErrL = calcrms(phaseNoise-phi)^2

bert = T_BER_SNR_mQAM(idbw(snr),mn);
h1=figure; plot(snr,bert,snr,ber); grid on
xlabel('SNR'); ylabel('BER');
h2=figure; plot(1:symlen,phaseNoise,'LineWidth',2); hold on
           plot(1:symlen,phi); grid on

mngFigureWindow(h1,h2);



% EOF
