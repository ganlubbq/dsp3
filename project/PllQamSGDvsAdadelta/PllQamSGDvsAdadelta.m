% COMPARE CPR USING PHASE-LOCK LOOP BASED ON SGD ALGORITHM AND ADADELTA
% ALGORITHM
%
% Decision directed PLL
%
% Framed data structure with training head
%
% Test cases:
%   see go.m
%
% M-QAM WITH TIME VARYING PHASE ERROR
function [ber] = PllQamSGDvsAdadelta(bitpersym, snr, pnvar, cfo, mu1, mu2)

if nargin < 1
    bitpersym = 2;
end
if nargin < 2
    snr = 20;
end
if nargin < 3
    pnvar = 1e-5;
end
if nargin < 4
    cfo = 1e-3;
end
if nargin < 5
    mu1 = 0.01;
end
if nargin < 6
    mu2 = 0.001;
end

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',283));

symlen = 2^18;
tvec = 0 : 1 : (symlen-1);
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
txLaserPnVar = 2 * pi * pnvar;
phaseNoise = genLaserPhaseNoise(symlen, txLaserPnVar, pi/6);
phaseNoise = phaseNoise(:);

% add phase noise with cfo
symTxPn = symTx .* exp(1i*phaseNoise) .* exp(1i*2*pi*cfo*tvec);

% add white gaussian noise
symTxPn = symTxPn + z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the nonlinear least squares problem using stochastic gradient
% descent algorithm. The cost function is J = |rx*exp(-i*phi)-tx|^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize stochastic gradient descent algorithm, implemented as a
% phase-lock loop, using 1 sample per symbol
symRec = zeros(size(symTxPn));
phi = zeros(size(symTxPn));
nco = zeros(size(symTxPn));
grad = zeros(size(symTxPn));
J = zeros(size(symTxPn));
symRec(1) = symTxPn(1);
trainingLength = 8;
frameLength = 512;
for k = 2 : length(symTxPn)
    % output
    symRec(k) = symTxPn(k) .* exp(-1i * phi(k-1));
    
    if mod(k, frameLength) <= trainingLength
        grad(k) = -imag(symRec(k) .* conj(symTx(k)));
    else
        grad(k) = -imag(symRec(k) .* conj(makeHardDecision(symRec(k), 2^bitpersym)));
    end
    
    % err integration
    nco(k) = nco(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient
    phi(k) = phi(k-1) - mu1 * grad(k) - mu2 * nco(k);
    
    % squared error
    J(k) = abs(symRec(k) - symTx(k)) .^ 2;
end

truePhase = unwrap(angle((symTxPn - z) .* conj(symTx)));

bitRx = slicerGrayQam(symRec, 2^bitpersym);
ber(1) = sum(abs(bitTx(:) - bitRx(:))) / (bitpersym * symlen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the nonlinear least squares problem using stochastic gradient
% descent algorithm with adaptive step size - Adadelta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symRec = zeros(size(symTxPn));
phi = zeros(size(symTxPn));
nco = zeros(size(symTxPn));
grad = zeros(size(symTxPn));
J = zeros(size(symTxPn));
symRec(1) = symTxPn(1);
gamma = 0.95;
meanSquareGradient = 0;
meanSquareDeltaPhi = 0;
epsilon = 1e-5;
for k = 2 : length(symTxPn)
    % output
    symRec(k) = symTxPn(k) .* exp(-1i * phi(k-1));
    
    if mod(k, frameLength) <= trainingLength
        grad(k) = -imag(symRec(k) .* conj(symTx(k)));
    else
        grad(k) = -imag(symRec(k) .* conj(makeHardDecision(symRec(k), 2^bitpersym)));
    end
    
    % update mean square of gradient
    meanSquareGradient = gamma * meanSquareGradient + (1 - gamma) * grad(k)^2;
    
    % rms of gradient
    rmsGradient = sqrt(meanSquareGradient + epsilon);
    
    % rms of delta phi
    rmsDeltaPhi = sqrt(meanSquareDeltaPhi + epsilon);
    
    % error
    deltaPhi(k) = (rmsDeltaPhi / rmsGradient) * grad(k);
    
    % update mean square of delta phi
    meanSquareDeltaPhi = gamma * meanSquareDeltaPhi + (1 - gamma) * deltaPhi(k)^2;
    
    % error integration
    nco(k) = gamma * nco(k-1) + (1 - gamma) * deltaPhi(k);
    
    % update filter coeff. along opposite direction of gradient
    phi(k) = phi(k-1) - deltaPhi(k);% - nco(k);
    
    % squared error
    J(k) = abs(symRec(k) - symTx(k)) .^ 2;
end

truePhase = unwrap(angle((symTxPn-z) .* conj(symTx)));

bitRx = slicerGrayQam(symRec, 2^bitpersym);
ber(2) = sum(abs(bitTx(:) - bitRx(:))) / (bitpersym * symlen);

return
