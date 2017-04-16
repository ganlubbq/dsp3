function H = calcBesselResponse(nSample, fs, order, bandwidth)
% Calculate frequency response of Bessel filter with input order and
% 3dB bandwidth
% 
% Example: 
%   H = calcBesselRespon(nSample, fs, order, bandwidth)
% 
% Input: 
%       nSample     - number of samples
%       fs          - sampling frequency [Hz]
%       order       - filter order
%       bandwidth   - filter 3dB bandwidth in [Hz]
% 
% Reference: 
% 
% Note: the intrinsic group delay is removed in this implementation
% 
% See also: calcDispResponse, calcOptGaussFlt
vAlpha = [1.0, 1.361654129, 1.755672389, 2.113917675, 2.427410702, 2.703395061];

% frequency grid in Hz
freqGrid = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] * fs / nSample;

% Laplace transformation
s = 1i * freqGrid / bandwidth * vAlpha(order);

% Order coeff
n = order;
for k = 0 : n
    b(k+1) = factorial(2*n-k) / (2^(n-k) * factorial(k) * factorial(n-k));
    D(k+1,:) = b(k+1) * s.^k;
end

% numerator/denominator
H = D(1,:) ./ sum(D);
H = H.';

% Remove group delay
ndxFreq = find(freqGrid >= 0 & freqGrid < bandwidth / vAlpha(order));
pp = polyfit(freqGrid(ndxFreq) / 1e9, unwrap(angle(H(ndxFreq))), 1);
H = H .* exp(-1i * polyval(pp, freqGrid/1e9));

% normalize
H = H / abs(H(1));

return
