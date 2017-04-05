% TEST SCRIPT FOR IMPLEMENT A 2ND-ORDER PHASE-LOCK LOOP BASED ON SGD ALGORITHM
%
% Note that one sample per symbol is used for the loop. In the case of
% large carrier frequency offset, it is necessary to turn on the 2nd-order
% loop filter. One can turn on and off the 2nd order loop filter to observe
% the phase tracking results.
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
% M-QAM WITH TIME VARYING PHASE ERROR
function [] = PllQamSGD(bitpersym, snr, mu1, mu2)

if nargin < 1
    bitpersym = 2;
end
if nargin < 2
    snr = 20;
end
if nargin < 3
    mu1 = 0.01;
end
if nargin < 4
    mu2 = 0.001;
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
% mu1 = 0.01;         % gain parameter 1st-order
% mu2 = 0.001 * 1;    % gain parameter 2nd-order
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
    phi(k) = phi(k-1) - mu1 * grad(k) - mu2 * nco(k);
    
    % squared error
    J(k) = abs(symRec(k) - symTx(k)) .^ 2;
end

truePhase = unwrap(angle((symTxPn-z) .* conj(symTx)));

% phi = mod(phi, 2*pi);
% truePhase = mod(truePhase, 2*pi);

figure; 
subplot(221); plot(symTxPn,'.'); grid on; %axis([-2.5 2.5 -2.5 2.5]);
subplot(222); plot(symRec,'.'); grid on; %axis([-2.5 2.5 -2.5 2.5]);
subplot(223); plot(tvec,phi,tvec,truePhase,'r'); grid on
subplot(224); plot(dbw(J(1:1000))); grid on; %xlim([0 symlen]); ylim([-100 20])
return
