function [y, filter_coef] = least_squares_filter(x, ts, alg, gain, lambda, filter_len)
% Linear adaptive filter with least squares principle.
% Two specific algorithms are implemented, i.e., least mean squares (LMS),
% recursive least squares (RLS).
%
% y = least_squares_filter(x, ts, alg, gain, lambda, filter_len)
%
% Inputs:
%   x : measured data to be filtered
%   ts : training sequence with same size of input 'x'
%   alg : name of algorithm, e.g. 'LMS', 'RLS'
%   gain : could be [] in RLS case
%   lambda : forgetting factor, could be [] in LMS case
%   filter_len : filter length
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
            filter_coef(:, ii + 1) = filter_coef(:, ii)...
                - gain * (xx.' * filter_coef(:, ii) - ts(ii)) ...
                * conj(xx) ./ (abs(xx) + eps);
            %%% or...
%             filter_coef(:, ii + 1) = filter_coef(:, ii)...
%                 - gain * (xx' * filter_coef(:, ii) - conj(ts(ii))) ...
%                 * xx ./ (abs(xx) + eps);
            y(ii) = filter_coef(:, ii + 1).' * xx;
        end
    case 'rls'
        for ii = 1 : length(x) - filter_len + 1
            xx = x_ext((1 : filter_len) + (ii - 1));
            gain = (1/lambda) * sigma * xx / (1.0 + (1/lambda) * xx' * sigma * xx);
            filter_coef(:, ii + 1) = filter_coef(:, ii)...
                + conj(gain) * (ts(ii) - xx.' * filter_coef(:, ii));
            sigma = (1/lambda) * (eye(filter_len) - gain * xx') * sigma;
            y(ii) = filter_coef(:, ii + 1).' * xx;
        end
    otherwise
        keyboard; 
end
return
