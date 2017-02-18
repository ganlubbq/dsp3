% DESCRIPTION
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
% Copyright 2015 default
function p = find_period(rx)

[c,lags] = xcorr(rx(1:1000),rx);

c = abs(c);

first_half = c(1001:end/2);
second_half = c(end/2+1000:end);

m1 = mean(first_half);
m2 = mean(second_half);

figure; plot(lags,c);

if m1>m2
    cc = first_half - m1;
else
    cc = second_half - m2;
end

idx = find(cc>max(cc)/2);

p = diff(idx);

return
