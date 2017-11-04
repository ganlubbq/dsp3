function [new_x] = zero_pad(x, zp_length, mode)
% Zero-pad the input data
%
% Call: [new_x] = zero_pad(x, zp_length, mode)
%
% Input:
%       x: column vector, input data
%       zp_length: scalar, zero-padding length
%       mode: string, 'front' or 'back'

if ~iscolumn(x), warning('input::column vector'); keyboard; end
if strcmpi(mode, 'front')
    new_x = [zeros(zp_length, size(x,2)); x];
elseif strcmpi(mode, 'back')
    new_x = [x; zeros(zp_length, size(x,2))];
else
    keyboard;
end
