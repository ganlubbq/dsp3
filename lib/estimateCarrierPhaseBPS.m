function [theta] = estimateCarrierPhaseBPS(observations, blocksize, M, mn)
% Blind phase search estimation of the carrier phase noise
% 
% Note:
%
% Example:
%   
% See also:

% How to solve the cycle-slip problem
phi = (0:M-1) / M * pi/2;

observations = normalizeQam(observations, mn);

% data matrix
xp = repmat(observations(:),1,M);

% test phase matrix
phsmat = ones(length(observations),1) * exp(1j*phi);

% rotated data
xpr = xp .* phsmat;

% distance between rotated data and its hard decision
for ii = 1:M
	D(:,ii) = smooth(abs(xpr(:,ii) - makeHardDecision(xpr(:,ii),mn)).^2, blocksize);
end

[~,idx] = min(D.');

% unwrap the phase estimation
theta = unwrap(4 * phi(idx)) / 4;

theta = theta(:);

return