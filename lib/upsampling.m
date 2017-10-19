function output = upsampling(input, sps, insert_zero)
% Upsampling signal by holding the sample value by default
% upsampling(input, sps, insert_zero)
% set 'insert_zero' to 1 to insert zeros between samples
if ~iscolumn(input)
    error('the first input has to be a column vector');
end
temp = ones(sps, 1) * input.';
if insert_zero
    temp(2 : end, :) = 0; % insert zeros
end
output = temp(:);
return
