% Test script for calculating theoretical coding rate of bianry symmetric
% channel with constrained input power, i.e., the channel capacity is
% related to the snr. The minimum required ebno for binary data transmission
% is found to be -1.6dB.
clear

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1));

% simulating a binary symmetric channel with a bpsk modulation
bitpersym = 1;
mn = 2^bitpersym;
symlen = 2^15;
a = constellation(mn);

% signal power
sp = sum(abs(a).^2) / mn;

% SNR with dB unit
snr = -20 : 0.5 : 10;

for pp = 1 : length(snr)
    % noise power with dB unit 
	th2db = 10*log10(sp) - snr(pp); 
	th2ln = 10^(th2db/10);
    
    % add gaussian real noise
	z = gaussian_noise(1, symlen, th2db, 'dbw', 'real');
	
	for kk = 1 : mn
		for ii = 1 : mn
			tmp(ii, :) = exp(-(abs(a(kk) + z - a(ii)).^2 - abs(z).^2) / (2*th2ln));
		end
		ss(kk, :) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss, 2));
	
    % simulated constrained capacity
	chc(pp) = log2(mn) - c2 / mn; 
end

% Limit for real modulation such as bpsk
lcc = 0.5 * log2(1 + idbw(snr));

% plot
h1 = figure(1); hold on; 
plot(snr, lcc, 'k', 'LineWidth', 2); 
plot(snr, chc, 'r', 'LineWidth', 2); 
grid on
ylim([0,2]);
xlabel('SNR (dB)'); 
ylabel('Capacity (bit/symbol)'); 
legend('Gaussian', 'BPSK');

% esno = .5*snr for real modulation
esno = 10*log10(0.5.*10.^(snr/10)); 

% theoretical code rate is equal to capacity
ebno = esno - 10*log10(bitpersym) - 10*log10(chc); 
h2 = figure; 
plot(ebno, chc, 'LineWidth', 2);
grid on; 
xlabel('E_b/N_0'); 
ylabel('Coding Rate');
legend('R < C = 1 - H(p)');

mngFigureWindow(h1,h2);

% EOF