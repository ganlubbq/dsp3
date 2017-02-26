function snr = osnr2snr(osnr, rs, sps, dmod)
%OSNR2SNR Convert OSNR in 0.1nm (12.5GHz) in dB to SNR in dB
% 
% Example: snr = osnr2snr(osnr,rs,sps,dmod)
% 
% Input: 
%       osnr     - osnr in db
%       rs       - symbol rate
%       sps      - sample per symbol
%       dmod     - data mode, complex or real

switch lower(dmod)
    case 'complex'
        snr = osnr - 10*log10(sps) - 10*log10(rs/12.5e9);
    case 'real'
        snr = osnr - 10*log10(sps) - 10*log10(rs/12.5e9) + 10*log10(2);
    otherwise
        warning('date mode not supported');
        keyboard;
end

return
