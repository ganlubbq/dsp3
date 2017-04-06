clear
snr = 5 : 0.5 : 10;
for ii = 1 : length(snr)
    ber = PllQamSGDvsAdadelta(2, snr(ii), (1/28)*1e-3, 0.0, 0.1, 0.0)
    bert(ii, :) = ber;
end

figure;
semilogy(bert);
grid on