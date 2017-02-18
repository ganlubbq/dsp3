% Calculate the root mean squared value of input
% 
% Example: Y = CALCRMS(X)
% 
% Input: X has to be a 1d vector, either a column or a row vector
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright: default

function Y = calcrms(X)

if ~iscolumn(X) && ~isrow(X)
    warning('input vector has to be a 1d vector'); keyboard;
end

Y = sqrt(sum(abs(X).^2)/length(X));

return


