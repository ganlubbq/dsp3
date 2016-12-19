% ideal timing recovery of the time shifted version of input data
%
% Example: output = clockRecIdeal(after,before)
% 
% Input: 
%       before  - act as the reference
%       after   - act as the time shifted version
% 
% Referene: I. Fatadin, D. Ives, and S. J. Savory, ¡°Compensation of Frequency
% Offset for Differentially Encoded 16- and 64-QAM in the Presence of Laser
% Phase Noise,¡± IEEE Photonics Technology Letters, vol. 22, no. 3, pp.
% 176¨C178, Feb. 2010.
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function output = clockRecIdeal(after,before)

% the longer the better
MAX_XCOR_LEN = 4096;

if length(before) < MAX_XCOR_LEN
    MAX_XCOR_LEN = length(before);
end

[vals,lags] = xcorr(after(1:MAX_XCOR_LEN),before(1:MAX_XCOR_LEN));

[~,ndx]     = max(vals);

output      = circshift(after(:),-lags(ndx));

output      = output.';

return


