% Calculate theoretical ber of Gray-coded square mqam in awgn channel
% 
% Input: 
%   snr : snr in linear scale for one dimension; 
%   k : bit per symbol;
% 
% See Also: ber2snr
function ber = snr2ber(snr, k)
M = 2.^k;
scalar = sqrt(1/(2/3*(M-1)));

if k == 1 % bpsk
	esno = snr / 2;
else
	esno = snr;
end

ser = 2 * (1-1/sqrt(M)) * erfc(scalar * sqrt(esno)) - ...
    (1-2/sqrt(M)+1/M) * (erfc(scalar * sqrt(esno))).^2;
             
if k == 1
    ber = 0.5 * erfc(sqrt(esno));
elseif k == 2
    ber = 1 - sqrt(1 - ser);
else % todo
    ber = ser / k;
end

return
