function estimated = estimatorLeastSquaresBatch(signals, observations)
%ESTIMATORLEASTSQUARESBATCH Implementing the least squares estimator based
%on a linear model.
% 
% Linear model: signals * parameter_to_be_estimated = observations
%
%       H * \theta = x
%
% Signals in the input is denoted by H and is a Nxp matrix. Observations is
% denoted by x and is a Nx1 vector. The estimation of model parameters is
% given by (H^T H)^-1 * H^T x and is a px1 vector.
%
% Generally, N is required to be greater than or equal to p to form an
% overdetermined system described by N linear equations with p unknowns.

if size(signals,1) < size(signals,2)
    warning('underdetermined system');
end

% solving an overdetermined system to minimize squared error
% estimated = inv(signals'*signals)*signals'*observations;
estimated = (signals'*signals)\(signals')*observations;

end

