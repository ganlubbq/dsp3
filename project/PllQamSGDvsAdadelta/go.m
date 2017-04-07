clear

snr = 5 : 1 : 14;

for ii = 1 : length(snr)
    ber = PllQamSGDvsAdadelta(2, snr(ii), (0.5/28)*1e-3, 2e-3, 0.02, 0.002)
    bert(ii, :) = ber;
end

figure;
semilogy(snr, bert, 'LineWidth', 2);
grid on