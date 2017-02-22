% Upsampling signal by inserting zeros
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright DEFAULT

function output = upSampInsertZeros(input, sps)

if ~iscolumn(input)
    error('the first input has to be a column vector');
end

temp = ones(sps,1) * input.';
temp(2:end,:) = 0; % insert zeros

output = temp(:);

return
