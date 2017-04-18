function H = calcRCFreqResponse(nSample, fs, fbaud, alpha, mode)
% CALCULATE RAISED COSINE DIGITAL FILTER RESPONSE IN FREQUENCY DOMAIN
% 
if nargin < 5
    mode = 'rc';
end

freq = [0:nSample/2]' / (nSample/2) * fs * 0.5;

H = zeros(size(freq));
f_low = (1-alpha) * fbaud/2;
f_high = (1+alpha) * fbaud/2;

fndxl = find(freq <= f_low);
fndxm = find(freq > f_low & freq <= f_high);

H(fndxl) = ones(size(fndxl));
H(fndxm) = 0.5 * (1 + cos(pi / alpha * (freq(fndxm) / fbaud - (1 - alpha) / 2)));

switch lower(mode)
    case 'rc'
        % rcos
    case 'rrc'
        H = sqrt(H);
    otherwise
        keyboard;
end

H = [H; conj(flipud(H(2:end-1)))];

return
