function esno = osnr2esno(osnr,rs)
%OSNR2ESNO Convert OSNR in 0.1nm (12.5GHz) in dB to ESNO in dB
% 

esno = osnr - 10*log10(rs/12.5e9);

return