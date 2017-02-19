%% TEST SCRIPT FOR CARRIER PHASE ESTIMATION WITH LMS FILTER OR PHASE-LOCK LOOP OF 1ST ORDER
% THE LEAST SQUARES CRITERIA IS ASSUMED. 
%
% Note that one sample per symbol is used for the loop. In the case
% of large carrier frequency offset, it is necessary to turn on the
% 2nd-order loop filter. One can turn on and off the 2nd order loop filter
% to observe the phase tracking results.
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
%% QPSK WITH TIME VARYING PHASE ERROR
clear

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

bitpersym = 2;

symlen = 2^16;
tvec = (0:(1/2e9):(symlen-1)) * (1/2e9);

% input power
a = constellation(2^bitpersym);
sp = sum(abs(a).^2) / (2^bitpersym);

% _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
snr = 13;
sigma2 = 10*log10(sp) - snr;
z = wgn(1,symlen,sigma2,'dbw','complex');

% modulation
bitTx = randi([0 1],bitpersym,symlen);
symTx = symbolizerGrayQam(bitTx);

% define phase noies
txLaserPnVar = 2*pi*10e3/2e9;
phaseNoise = genLaserPhaseNoise(symlen,txLaserPnVar,pi/6);

% add phase noise with cfo
cfo = 20e6;
symTxPn = symTx .* exp(1j*phaseNoise) .* exp(1j*2*pi*cfo*tvec);

% add white gaussian noise
symTxPn = symTxPn + z;


%% Solving the nonlinear LS function with gradient descent method
% Using least squares equalization model, i.e., J = |x*exp(-j*phi)-ref|^2

% initialize stochastic gradient descent algorithm, implemented as a
% phase-lock loop, using 1 sample per symbol
mu1 = 0.09;  % gain parameter 1st-order
mu2 = 0.001 * 1; % gain parameter 2nd-order
symRec(1) = symTxPn(1);
phi(1) = 0;
nco(1) = 0;
for k = 2:length(symTxPn)
    % output
    symRec(k) = symTxPn(k) .* exp(-1j * phi(k-1));
    % stochastic gradient, also a PED with a typical S-curve
    grad(k) = -imag(symRec(k) .* conj(symTx(k)));
     % err integration
    nco(k) = nco(k-1) + grad(k);
    % update filter coeff. along opposite direction of gradient
    phi(k) = phi(k-1) - mu1*grad(k) - mu2*nco(k);
    % squared error
    J(k) = abs(symRec(k) - symTx(k)).^2;
end

truePhase = unwrap(angle((symTxPn-z) .* conj(symTx)));

figure; 
subplot(221); plot(symTxPn,'.'); grid on; axis([-2.5 2.5 -2.5 2.5]);
subplot(222); plot(symRec,'.'); grid on; axis([-2.5 2.5 -2.5 2.5]);
subplot(223); plot(tvec,phi,tvec,truePhase,'r'); grid on
subplot(224); plot(dbw(J)); grid on; xlim([0 symlen]); ylim([-100 20])



% EOF
