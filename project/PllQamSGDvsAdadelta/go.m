clear

for snr = 11 : 20
    ber = PllQamSGDvsAdadelta(2, snr, 2e-4, 0.0, 0.01, 0.0)
    bert(snr, :) = ber;
end

figure;
semilogy(bert);
grid on