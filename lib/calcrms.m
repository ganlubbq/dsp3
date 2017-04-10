function Y = calcrms(X)
% Calculate the root mean squared value of input
% 
% Example: Y = CALCRMS(X)
% 
% Input: X has to be a 1d vector, either a column or a row vector
% 
% See Also: 
% 
if ~iscolumn(X) && ~isrow(X)
    warning('input vector has to be a 1d vector'); keyboard;
end

Y = sqrt(sum(abs(X).^2) / numel(X));

return
