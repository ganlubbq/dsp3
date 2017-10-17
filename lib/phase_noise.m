function phase_noise = genLaserPhaseNoise(nSample,pnvar,p0)
% Generate laser phase noise with certain variance
% 
% Example: phase_noise = genLaserPhaseNoise(nSample,pnvar,p0)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 DEFAULT

tmp = randn(1,nSample);

% remove dc component
tmp = tmp - mean(tmp);

% normalize
tmp = tmp ./ calcrms(tmp);

phase_noise = p0 + cumsum(tmp .* sqrt(pnvar));

return
