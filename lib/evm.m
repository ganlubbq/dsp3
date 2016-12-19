% Get error vector magnitude for input signals
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

function v = evm(x,mn)

xn = normalize_mqam(x,mn);
xd = graySlicer_mQAM(xn,mn);

Ps = mean(abs(xd).^2);
c = constellation(mn);

for ii = 1:mn
    idx{ii} = find(xd==c(ii));
    n(ii) = mean(abs(xn(idx{ii})-c(ii)).^2);
end

v = mean(n)/Ps;
% v = sqrt(v);
return

% function c = constellation(mn)
% switch mn
%     case 1
%         c = [complex(-0,0); complex(1,0)];
%     case 2
%         c = [complex(-1,0); complex(1,0)];
%     case 4
%         x = [+1 -1]; quar = [x;x]; quar = quar(:);
%         y = [-1 +1]; inph = [y y]; inph = inph.';
%         c = inph + 1j*quar;
%     case 16
%         x = [3 1 -1 -3]; quar = [x;x;x;x]; quar = quar(:);
%         y = [-3 -1 1 3]; inph = [y y y y]; inph = inph.';
%         c = inph + 1j*quar;
%     case 64
%         x = [7 5 3 1 -1 -3 -5 -7]; quar = [x;x;x;x;x;x;x;x]; quar = quar(:);
%         y = [-7 -5 -3 -1 1 3 5 7]; inph = [y y y y y y y y]; inph = inph.';
%         c = inph + 1j*quar;
% end

