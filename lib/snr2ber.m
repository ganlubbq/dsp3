% Calculate theoretical ber of Gray-coded square mqam in awgn channel. The
% complex modulated square mqam signal can be considered as a 2-dimensional
% signal and the symbol error rate can be expressed as 
% P_sym = 1 - (1 - P_i/q)^2 = 2 * P_i/q - (P_i/q)^2,
% where P_i/q is the symbol error rate of each quadrature
% 
% Input: 
%   snr : snr in one dimension; 
%   k : bit-per-symbol;
% 
% See Also: ber2snr
function ber = snr2ber(snr, k, snrmode)
if nargin < 3
    warning('You have to specify the SNR unit (dB or linear)');
    keyboard;
end

M = 2.^k;
scalar = sqrt(1/(2/3*(M-1)));

switch lower(snrmode)
    case 'db'
        snr = 10 .^ (0.1 * snr);
    case 'linear'
        % 
    otherwise
        keyboard;
end

if k == 1 % bpsk
	esno = snr / 2;
else
	esno = snr;
end

ser_iq = (1 - 1/sqrt(M)) * erfc(scalar * sqrt(esno));

ser = 2 * ser_iq - ser_iq.^2;

if k == 1
    ber = 0.5 * erfc(sqrt(esno));
elseif k == 2
    ber = 1 - sqrt(1 - ser);
else % todo
    ber = ser / k;
end

return
