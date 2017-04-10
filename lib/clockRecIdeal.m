function output = clockRecIdeal(after,before)
% ideal timing recovery of the time shifted version of input data
%
% Example: output = clockRecIdeal(after,before)
% 
% Input: 
%       before  - act as the reference
%       after   - act as the time shifted version
% 
% Note: 
% 
% See Also: 
% the longer the better
MAX_XCOR_LEN = 4096;

if length(before) < MAX_XCOR_LEN
    MAX_XCOR_LEN = length(before);
end

[vals,lags] = xcorr(after(1:MAX_XCOR_LEN), before(1:MAX_XCOR_LEN));

[~,ndx]     = max(vals);

output      = circshift(after(:), -lags(ndx));

output      = output.';

return
