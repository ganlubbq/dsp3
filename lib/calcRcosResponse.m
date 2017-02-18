% CALCULATE RAISED COSINE DIGITAL FILTER RESPONSE IN FREQUENCY DOMAIN
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
% Copyright Default

function H = calcRcosResponse(nSample,fs,fbaud,alpha,mode)

if nargin<5
    mode = 0;
end

freq = [0:nSample/2]'/(nSample/2)*fs*0.5;

H = zeros(size(freq));
f_low = (1-alpha)*fbaud/2;
f_high = (1+alpha)*fbaud/2;

fndxl = find(freq<=f_low);
fndxm = find(freq>f_low & freq<=f_high);
fndxh = [fndxl;fndxm];

H(fndxl) = ones(size(fndxl));
H(fndxm) = 1/2*(1+cos(pi/alpha*(freq(fndxm)/fbaud-(1-alpha)/2)));

% for ii = 1:length(f)
%     if abs(f(ii)) <= (1-beta)/(2*T)
%         H(ii) = T;
%     elseif (abs(f(ii))>(1-beta)/(2*T) && abs(f(ii))<=(1+beta)/(2*T))
%         H(ii) = (T/2) * (1+cos(pi*T/beta*(abs(f(ii))-(1-beta)/(2*T))));
%     else
%         H(ii) = 0;
%     end
% end

switch mode
    case 0
        % rcos, do null
    case 1
        H = sqrt(H);
    otherwise
        error('raised cosine configuration not supported');
end

H = [H;conj(flipud(H(2:end-1)))];

return


