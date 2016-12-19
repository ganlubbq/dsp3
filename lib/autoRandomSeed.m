% Reset system random number seed automatically
% 
% Example: autoRandomSeed(seed)
% 
% Input: rnsd - seed
% 
% Reference: 
% 
% Note: 
% 
% See Also: rand, randi
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 
function [] = autoRandomSeed(rnsd)

if rnsd
    rng(rnsd);
else
    rng_now = rng;
    tmp = randi([1,65534],1,1);
    while (tmp == rng_now.Seed)
        tmp = randi([1,65534],1,1);
    end
    rng(tmp);
end

return
