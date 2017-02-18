%T_BER Gray coded square M-QAM theoretical ber with symbol rate of Rs
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
% Copyright default

function t_ber = T_BER_mQAM(dB_osnr,M,Rs)

scale_factor    = sqrt(1/(2/3*(M-1)));
k               = log2(M);

OSNR            = 10.^(dB_osnr/10);

EsNo            = OSNR * 12.5e9 / Rs;
EsNo_db         = 10 * log10(EsNo);

EbNo            = EsNo / k;
EbNo_db         = 10 * log10(EbNo);

t_ser           = 2 * (1-1/sqrt(M)) * erfc( scale_factor * sqrt(EsNo) ) - ...
                 (1-2/sqrt(M)+1/M) * (erfc(scale_factor * sqrt(EsNo))).^2;
if M==2
    t_ber = 0.5*erfc(sqrt(EsNo));
elseif M==4
    t_ber = 1-sqrt(1-t_ser);
else
    t_ber           = t_ser / k;
end

return


