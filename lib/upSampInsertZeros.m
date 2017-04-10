function output = upSampInsertZeros(input, sps)
% Upsampling signal by inserting zeros
% 
% See Also: upSampInsertZeros
if ~iscolumn(input)
    error('the first input has to be a column vector');
end

temp = ones(sps,1) * input.';
temp(2:end,:) = 0; % insert zeros

output = temp(:);

return
