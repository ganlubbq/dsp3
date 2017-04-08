function theta = estimateCarrierPhaseML(observations, blocksize, mn)
% Implementing the carrier phase estimation based on the decision-directed
% maximum likelihood criterial. 
%
% The algorithm works blockwise, within each block the phase is assumed to
% be deterministic and constant.
%
% mn in the input stands for mn-QAM
if ~iscolumn(observations)
    error('first input has to be a column vector');
end

% based on decision
signals = makeHardDecision(observations, mn);

% sliding window smoothing
H = filter(ones(1,blocksize)/blocksize, 1, observations .* conj(signals));

% the use of angle func. limits the estimation range
theta = angle(H);

return
