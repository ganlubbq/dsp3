function [y, filter_coef] = least_squares_filter(x, ts, mf, alg, gain, lambda, filter_len)
% Linear adaptive filter with least squares principle.
% Two specific algorithms are implemented, i.e., least mean squares (LMS),
% recursive least squares (RLS).
%
% y = least_squares_filter(x, ts, mf, alg, gain, lambda, filter_len)
%
% Inputs:
%   x : column vector, measured data to be filtered 
%   ts : column vector, training sequence with same size of input 'x', 
%        when tx is NAN, the filter will switch to decision-directed mode
%   mf : string, modulation format, e.g. 'QPSK', '16QAM'
%   alg : string, name of algorithm, e.g. 'LMS', 'RLS'
%   gain : scalar, gain factor, could be [] in RLS case
%   lambda : scalar, forgetting factor, could be [] in LMS case
%   filter_len : scalar, filter-tap length
%
% Note the numerical stability problem of conventional RLS algorithm
% Note the normalization realization of LMS algorithm

%%% check the size of inputs
if length(x) ~= length(ts), keyboard; end
if ~iscolumn(x), warning('input::column vector'); keyboard; end

%%% extend measured data
x_ext = [zeros(filter_len - 1, 1); x];

%%% initialization - filter coefficients
filter_coef = zeros(filter_len, length(x) - filter_len + 2);
%%% initialization - covariance matrix
sigma = 1e5 * eye(filter_len);

y = zeros(size(x));
switch lower(alg)
    case 'lms'
        for ii = 1 : length(x) - filter_len + 1
            xx = x_ext((1 : filter_len) + (ii - 1));
            y(ii) = xx.' * filter_coef(:, ii);
            if isnan(ts(ii))
                target = hard_decision(y(ii), mf2mn(mf));
            else
                target = ts(ii);
            end
            filter_coef(:, ii + 1) = filter_coef(:, ii)...
                - gain * (y(ii) - target) ...
                * conj(xx) ./ (abs(xx) + eps);
            %%% alternative
%             y(ii) = xx' * filter_coef(:, ii);
%             if isnan(ts(ii))
%                 target = hard_decision(y(ii));
%             else
%                 target = ts(ii);
%             end
%             filter_coef(:, ii + 1) = filter_coef(:, ii)...
%                 - gain * (y(ii) - conj(target)) ...
%                 * xx ./ (abs(xx) + eps);
        end
    case 'rls'
        for ii = 1 : length(x) - filter_len + 1
            xx = x_ext((1 : filter_len) + (ii - 1));
            gain = (1/lambda) * sigma * xx / (1.0 + (1/lambda) * xx' * sigma * xx);
            y(ii) = xx.' * filter_coef(:, ii);
            if isnan(ts(ii))
                target = hard_decision(y(ii), mf2mn(mf));
            else
                target = ts(ii);
            end
            filter_coef(:, ii + 1) = filter_coef(:, ii)...
                + conj(gain) * (target - y(ii));
            sigma = (1/lambda) * (eye(filter_len) - gain * xx') * sigma;
        end
    otherwise
        keyboard; 
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
