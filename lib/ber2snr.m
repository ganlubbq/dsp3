% Calculate required snr for a given ber with mqam in the awgn model
%
% Algorithm:
%   qpsk : ber = 1 - sqrt(ser);
%          ser = erfc(sqrt(snr/2)) - (1/4)erfc^2(sqrt(snr/2));
%
% Note: erfc = 1 - erf
%       erfinv(erf(x)) = x
%
% Test: ber2snr(snr2ber(1.5,2,'db'),2,'db')
function snr = ber2snr(ber, k, snrmode)
if k == 1
    snr = 2 * (erfinv(1 - 2*ber)) .^ 2;
elseif k == 2
    snr = 2 * (erfinv(1 - 2*ber)) .^ 2;
else
    % todo: may use interpolation
end

switch lower(snrmode)
    case 'db'
        snr = 10 * log10(snr);
    case 'linear'
        % 
    otherwise
        keyboard;
end
return
