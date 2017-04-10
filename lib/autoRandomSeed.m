function [] = autoRandomSeed(rnsd)
% Reset system random number seed automatically
% 
% Example: autoRandomSeed(seed)
% 
% Input: rnsd - seed
% 
% Note: 
% 
% See Also: rand, randi
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
