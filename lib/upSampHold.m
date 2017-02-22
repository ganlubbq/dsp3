% Upsampling signal by holding the symbol value
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: upSampInsertZeros
% 
% Copyright DEFAULT

function output = upSampHold(input, sps)

if ~iscolumn(input)
    error('the first input has to be a column vector');
end

temp = ones(sps,1) * input.';

output = temp(:);

return
