function ber = T_BER_SNR_mQAM(SNR,M)
%T_BER Gray coded square M-QAM theoretical ber
% 
% Example: 
% 
% Input: SNR in linear scale for one dimension; M the alphabet size
% 
% Reference: 
% 
% Note: 
% 
% See Also: 

scale_factor    = sqrt(1/(2/3*(M-1)));
k               = log2(M);

if k == 1 % bpsk
	EsNo    = SNR /2;
else
	EsNo    = SNR;
end

EsNo_db         = 10 * log10(EsNo);
EbNo            = EsNo / k;
EbNo_db         = 10 * log10(EbNo);

t_ser           = 2 * (1-1/sqrt(M)) * erfc( scale_factor * sqrt(EsNo) ) - ...
                 (1-2/sqrt(M)+1/M) * (erfc(scale_factor * sqrt(EsNo))).^2;
             
if M == 2
    ber = 0.5*erfc(sqrt(EsNo));
elseif M==4
    ber = 1-sqrt(1-t_ser);
else % todo
    ber = t_ser / k;
end

return


