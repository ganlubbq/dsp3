% Calculate required snr for a given ber with mqam in the awgn model
%
% Algorithm:
%   qpsk : ber = 1 - sqrt(ser);
%          ser = erfc(sqrt(snr/2)) - (1/4)erfc^2(sqrt(snr/2));
%
% Note: erfc = 1 - erf
%       erfinv(erf(x)) = x
function snr = ber2snr(ber, k)
if k == 1
    snr = 2 * (erfinv(1 - 2*ber)) .^ 2;
elseif k == 2
    snr = 2 * (erfinv(1 - 2*ber)) .^ 2;
end
return
