function theta = estimateCarrierPhaseML(signals, observations, blocksize)
% Implementing the carrier phase estimation based on the maximum likelihood
% criterial. The algorithm works blockwise, within each block the phase is
% assumed to be deterministic and constant. This algorithm should be ideal
% (when the block size is set to be 1) for purely phase-modelated signal
% such as m-psk
if ~iscolumn(observations)
    error('first input has to be a column vector');
end

% sliding window smoothing
H = filter(ones(1,blocksize)/blocksize, 1, observations .* conj(signals));

% the use of angle func. limits the estimation range
theta = angle(H);

return
