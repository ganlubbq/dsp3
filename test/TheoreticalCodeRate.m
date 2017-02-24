%% TEST SCRIPT FOR CALCULATING THEORETICAL CODING RATE
% THE THEORETICAL CODING RATE IS ALWAYS LESS THAN THE CHANNEL CAPACITY,
% WHICH IS A BINARY REAL-VALUED CHANNEL
%
%% BPSK as binary system
clear
clc

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',199));
% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

bitpersym = 1;
mn = 2^bitpersym;

symlen = 2^13;

a = constellation(mn);

% input power
sp = sum(abs(a).^2) / mn;

snr = -20:1:11; % dB unit

for pp = 1:length(snr)

    % dB unit noise power
	th2db = 10*log10(sp) - snr(pp); 
	th2ln = 10^(th2db/10);
    
	z = genWGN(1,symlen,th2db,'dbw','real');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2 - abs(z).^2) / (2*th2ln));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
    % simulated constrained capacity for binary channel
	cc(pp) = log2(mn) - c2/mn; 
end

%% Limit for real modulation such as bpsk
lcc = 0.5*log2(1 + idbw(snr));

%% Figure
h1 = figure(1); hold on; grid on
plot(snr,cc,'r.-'); plot(snr,lcc,'k-'); ylim([0,2]);
xlabel('SNR (dB)'); ylabel('Capacity (bit/symbol)'); 
legend('BPSK','Shannon');

%% Coding rate
% esno = .5*snr for real modulation
esno = 10*log10(0.5.*10.^(snr/10)); 
% theoretical code rate is equal to capacity
ebno = esno - 10*log10(bitpersym) - 10*log10(cc); 
h2=figure; grid on; 
plot(ebno,cc);
xlabel('E_b/N_0'); ylabel('Coding Rate');

mngFigureWindow(h1,h2);

% EOF