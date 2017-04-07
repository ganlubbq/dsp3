clear

snr = 0 : 1 : 10;

for ii = 1 : length(snr)
    ber = PllQamSGDvsAdadelta(2, snr(ii), (1.5/28)*1e-3, 0.0, 0.02, 0.0)
    bert(ii, :) = ber;
end

figure;
semilogy(snr, bert, 'LineWidth', 2);
grid on