function theta = estimateCarrierPhaseML(signals, observations, blocksize)
%ESTIMATECARRIERPHASEML Implementing the carrier phase estimation based on
%the decision directed ML criterial. The algorithm works blockwise, within
%each block the phase is assumed to be deterministic and constant. 
%

if ~iscolumn(signals)
    error('first input has to be a column vector');
end
if ~iscolumn(observations)
    error('second input has to be a column vector');
end

% sliding window smoothing the 
H = smooth(observations .* conj(signals), blocksize);

% the use of angle func. limits the estimation range; could output H
% directly and use it as compensation
theta = angle(H);

return
