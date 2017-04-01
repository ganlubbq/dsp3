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

% m-QAM WITH TIME VARYING PHASE ERROR
function [] = PhaseLockLoopQamDemo(bitpersym, snr)

if nargin < 1
    bitpersym = 2;
end
if nargin < 2
    snr = 20;
end

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

symlen = 2^16;
tvec = 0 : (1/2e9) : ((symlen-1) * (1/2e9));
tvec = tvec(:);

% input power
a = constellation(2^bitpersym);
sp = sum(abs(a).^2) / (2^bitpersym);

% modulation
bitTx = randi([0 1],bitpersym,symlen);
symTx = symbolizerGrayQam(bitTx);

% _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
sigma2 = 10*log10(sp) - snr;
z = wgn(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');

% define phase noies
txLaserPnVar = 2 * pi * 10e3 / 2e9;
phaseNoise = genLaserPhaseNoise(symlen, txLaserPnVar, pi/6);
phaseNoise = phaseNoise(:);

% add phase noise with cfo
cfo = 20e6;
symTxPn = symTx .* exp(1i*phaseNoise) .* exp(1i*2*pi*cfo*tvec);

% add white gaussian noise
symTxPn = symTxPn + z;

% Solving the nonlinear LS function with gradient descent method
% Using least squares equalization model, i.e., J = |x*exp(-j*phi)-ref|^2

% initialize stochastic gradient descent algorithm, implemented as a
% phase-lock loop, using 1 sample per symbol
mu1 = 0.01;         % gain parameter 1st-order
mu2 = 0.0001 * 1;    % gain parameter 2nd-order
symRec(1) = symTxPn(1);
phi(1) = 0;
nco(1) = 0;
for k = 2:length(symTxPn)
    % output
    symRec(k) = symTxPn(k) .* exp(-1i * phi(k-1));
    
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

% phi = mod(phi, 2*pi);
% truePhase = mod(truePhase, 2*pi);

figure; 
subplot(221); plot(symTxPn,'.'); grid on; %axis([-2.5 2.5 -2.5 2.5]);
subplot(222); plot(symRec,'.'); grid on; %axis([-2.5 2.5 -2.5 2.5]);
subplot(223); plot(tvec,phi,tvec,truePhase,'r'); grid on
subplot(224); plot(dbw(J)); grid on; xlim([0 symlen]); ylim([-100 20])
return
