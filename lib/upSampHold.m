function output = upSampHold(input, sps)
% Upsampling signal by holding the symbol value
% 
% See Also: upSampInsertZeros
if ~iscolumn(input)
    error('the first input has to be a column vector');
end

temp = ones(sps,1) * input.';

output = temp(:);

return
