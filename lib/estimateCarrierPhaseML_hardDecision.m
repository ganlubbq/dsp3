function theta = estimateCarrierPhaseML_hardDecision(observations, blocksize, mn)
%ESTIMATECARRIERPHASEML Implementing the carrier phase estimation based on
%the decision directed ML criterial. The algorithm works blockwise, within
%each block the phase is assumed to be deterministic and constant. 
%
% mn in the input stands for mn-QAM



% based on decision
signals = makeHardDecision(observations, mn);

% sliding window smoothing
H = smooth(observations(:) .* conj(signals), blocksize);

% the use of angle func. limits the estimation range
theta = angle(H);

return


