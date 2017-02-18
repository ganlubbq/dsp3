% DELAYINTERFEROMETER Delay input by several samples and beat with itself
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
% Copyright 2015 Default

function xout = delayInterferSample(xin, sps)

xin = xin(:);

% delay and interferring
dx = circshift(xin, -sps);
dy = xin.*conj(dx);
% ry = real(dy);
ry = reshape(dy,sps,[],size(dy,2));
if sps > 1
    xout = squeeze( sum(ry) );
else
    xout = squeeze( ry );
end
% figure; plot(xout(:,1),'.'); grid on

return


