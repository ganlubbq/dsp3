function theta = estimateCarrierPhaseDAML(observations, blocksize, mn)
%ESTIMATECARRIERPHASEML Implementing the carrier phase estimation based on
%the decision directed ML criterial. The algorithm works blockwise, within
%each block the phase is assumed to be deterministic and constant. 
%
% mn in the input stands for mn-QAM

if ~iscolumn(observations)
    error('first input has to be a column vector');
end

% based on decision
signals = makeHardDecision(observations, mn);

% sliding window smoothing
H = smooth(observations .* conj(signals), blocksize);

% the use of angle func. limits the estimation range
theta = angle(H);

return
