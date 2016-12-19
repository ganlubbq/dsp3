% Single polarization equalization using lms algorithm with multiple
% iterations
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
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 
function [yout,mse,deth,h1] = lms_sng_eq( xin,mn,sps,mu,ntaps,iter,h1)

x = normalizeQam(xin, mn);

% scale x to max{1,1} form
x = x /(sqrt(mn)-1);

% make sure the tap number is odd
ntaps = ntaps + ~mod(ntaps,2);

% taps initialization
if nargin < 7
    h1 = zeros(ntaps,1);
    h1(floor(ntaps/2)+1) = 1;
end

if length(h1) ~= ntaps
    error('filter taps length error');
end

% reference constellation points
cstl = constellation(mn)/(sqrt(mn)-1);

% pad few samples at begining and end of x with circular periodic condition
halfnt = floor(ntaps/2);
extendx = [ x(end-halfnt+1:end,:); x; x(1:halfnt,:)];

for ii = 1:iter
    [xx,mse,deth] = LMS_FILTER_sng(extendx,h1,ntaps,mu,cstl,sps);
    %%[xx,mse,deth] = MCMA_FILTER_sng(extendx,h1,ntaps,mu,radius,sps,errid);
end

% format the output
yout    = xx(1:sps:end,:);
mse     = mse(1:sps:end);
deth    = deth(1:sps:end);

return
