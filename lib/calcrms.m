% Calculate the root mean squared value of input
% 
% Example: Y = CALCRMS(X)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function Y = calcrms(X)

Y = sqrt(sum(abs(X).^2)/length(X));

return


