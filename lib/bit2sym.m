function sym = bit2sym(bit,mn,format)
% BIT2SYM Convert binary bits to mQAM symbols
% 
% Example: sym = bit2sym(bit,mn,format)
% 
% Input: 
%       format - specify the leftmost bit being LSB or MSB
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright (c) 2015 Dawei Wang 

n = log2(mn);

if strcmp(format,'MSB')
    bit = fliplr(bit);
    for k = 1:n
        bit(:,k) = bit(:,k) * (2^(k-1));
    end
elseif strcmp(format,'LSB')
    for k = 1:n
        bit(:,k) = bit(:,k) * (2^(k-1));
    end
else
    error('input bit format is not supported!!!');
end

sym = sum(bit,2) + 1;

return


