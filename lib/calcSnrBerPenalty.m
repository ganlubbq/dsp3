function [ penalty] = calcSnrBerPenalty( SNR, BER1, BER2, tarBER)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% everything must be in column and log10
if ~iscolumn(SNR)
    warning('input snr has to be a column vector'); keyboard;
end
if ~iscolumn(BER1)
    warning('input ber1 has to be a column vector'); keyboard;
end
% Comment out to support matrix form of BER2
% if ~iscolumn(BER2)
%     warning('input ber2 has to be a column vector'); keyboard;
% end
if ~isscalar(tarBER)
    warning('input tarBER has to be a scalar'); keyboard;
end

% set the resolution to be 0,01 dB
snrFine = linspace(SNR(1),SNR(end),(SNR(end)-SNR(1))/0.01+1);

ber1Fine = interp1(SNR,BER1,snrFine(:));
ber2Fine = interp1(SNR,BER2,snrFine(:));

[temp, ndx1] = min(abs(ber1Fine-tarBER));
[temp, ndx2] = min(abs(ber2Fine-tarBER));

penalty = snrFine(ndx2) - snrFine(ndx1);

end

