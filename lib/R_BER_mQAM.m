%R_BER Required OSNR for Gray coded square M-QAM with ber of ber
% 
% Example: osnr = R_BER_mQAM(ber,M,Rs)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2011 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function dB_osnr = R_BER_mQAM(ber,M,Rs)

scale_factor    = sqrt(1/(2/3*(M-1)));
k               = log2(M);

osnr_v          = 10:0.0001:40;
OSNR            = 10.^(osnr_v/10);

EsNo            = OSNR * 12.5e9 / Rs;
EsNo_db         = 10 * log10(EsNo);

EbNo            = EsNo / k;
EbNo_db         = 10 * log10(EbNo);

t_ser           = 2 * (1-1/sqrt(M)) * erfc( scale_factor * sqrt(EsNo) ) - ...
                 (1-2/sqrt(M)+1/M) * (erfc(scale_factor * sqrt(EsNo))).^2;
             
t_ber           = t_ser / k;

[~,idx]         = min(abs(t_ber-ber));

dB_osnr         = osnr_v(idx);

return

