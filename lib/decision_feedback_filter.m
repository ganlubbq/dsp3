function [y, flt_coef] = decision_feedback_filter(x, ts, mf, alg, gain, lambda, ff_len, fb_len)
% Decision feedback equalizer (DFE) with least squares principle.
%
% Call: [y, flt_coef] = decision_feedback_filter(x, ts, mf, alg, gain, lambda, ff_len, fb_len)
%
% Inputs:
%   x : column vector, measured data to be filtered
%   ts : column vector, training sequence with same size of input 'x', 
%        when tx is NAN, the filter will switch to decision-directed mode
%   mf : string, modulation format, e.g. 'QPSK', '16QAM'
%   alg : string, name of algorithm, e.g. 'LMS', 'RLS'
%   gain : scalar, gain factor, could be [] in RLS case
%   lambda : scalar, forgetting factor, could be [] in LMS case
%   ff_len : scalar, feedforward filter-tap length
%   fb_len : scalar, feedback filter-tap length
%
% Note the numerical stability problem of conventional RLS algorithm
% Note the normalization realization of LMS algorithm

%%% check the size of inputs
if length(x) ~= length(ts), keyboard; end
if ~iscolumn(x), keyboard; end

%%% extend measured data
x_ext = zero_pad(x, filter_len-1, 'front');
%%% initialization - filter coefficients
ff_coef = zeros(ff_len, length(x) - filter_len + 2);

y = zeros(size(x));
for ii = 1 : length(x) - filter_len + 1
    
end
flt_coef = [];
end

%%%%%%%%%%%%
function mn = mf2mn(mf)
switch lower(mf)
    case 'bpsk'
        mn = 2;
    case 'qpsk'
        mn = 4;
    case '16qam'
        mn = 16;
    case '32qam'
        mn = 32;
    case '64qam'
        mn = 64;
    otherwise
        warning('filter::modulation format'); keyboard;
end
return

return
