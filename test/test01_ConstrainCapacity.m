%% TEST SCRIPT FOR CALCULATING CONSTRAINED CHANNEL CAPACITY
% THE CHANNEL CAPACITY FOR A DISCRETE MEMORYLESS GAUSSIAN CHANNEL
%
% Ungerboeck G. Channel coding with multilevel/phase signals[J]. 
% IEEE Transactions on Information Theory,
% 1982, 28(1):55-67.
%
%% QPSK

clear
clc
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',199));
% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
bitpersym = 2;
mn = 2^bitpersym;

symlen = 2^13;

a = constellation(mn);

% input power
sp = sum(abs(a).^2)/mn;

snr = 0:1:36;

for pp = 1:length(snr)

	th2 = 10*log10(sp) - snr(pp);
	
    % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
	z = wgn(1,symlen,th2,'dbw','complex');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2-abs(z).^2)/10^(th2/10));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
	cc(pp) = log2(mn)-c2/mn;
	
end

h1 = figure(1); plot(snr,cc,'k.-'); grid on; hold on;

%% 8-QAM

bitpersym = 3;
mn = 2^bitpersym;

a = constellation(mn);

% input power
sp = sum(abs(a).^2)/mn;

snr = 0:1:36;

for pp = 1:length(snr)

	th2 = 10*log10(sp) - snr(pp);
	
    % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
	z = wgn(1,symlen,th2,'dbw','complex');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2-abs(z).^2)/10^(th2/10));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
	cc(pp) = log2(mn)-c2/mn;
	
end

h1 = figure(1); plot(snr,cc,'c.-'); grid on; hold on;

%% 16-QAM
bitpersym = 4;
mn = 2^bitpersym;

a = constellation(mn);

% input power
sp = sum(abs(a).^2)/mn;

snr = 0:1:36;

for pp = 1:length(snr)

	th2 = 10*log10(sp) - snr(pp);
	
    % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
	z = wgn(1,symlen,th2,'dbw','complex');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2-abs(z).^2)/10^(th2/10));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
	cc(pp) = log2(mn)-c2/mn;
	
end

h1 = figure(1); plot(snr,cc,'g.-'); grid on;

%% 32-QAM
bitpersym = 5;
mn = 2^bitpersym;

a = constellation(mn);

% input power
sp = sum(abs(a).^2)/mn;

snr = 0:1:36;

for pp = 1:length(snr)

	th2 = 10*log10(sp) - snr(pp);
	
    % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
	z = wgn(1,symlen,th2,'dbw','complex');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2-abs(z).^2)/10^(th2/10));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
	cc(pp) = log2(mn)-c2/mn;
	
end

h1 = figure(1); plot(snr,cc,'b.-'); grid on;

%% 64-QAM
bitpersym = 6;
mn = 2^bitpersym;

a = constellation(mn);

% input power
sp = sum(abs(a).^2)/mn;

snr = 0:1:36;

for pp = 1:length(snr)

	th2 = 10*log10(sp) - snr(pp);
	
    % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
	z = wgn(1,symlen,th2,'dbw','complex');
	
	for kk = 1:mn
		for ii = 1:mn
			tmp(ii,:) = exp(-(abs(a(kk)+z-a(ii)).^2-abs(z).^2)/10^(th2/10));
		end
		ss(kk,:) = log2(sum(tmp));
	end
	
	c2 = sum(mean(ss,2));
	
	cc(pp) = log2(mn)-c2/mn;
	
end

%% Limit for complex modulation
lcc = log2(1+idbw(snr));

h1 = figure(1); plot(snr,cc,'r.-'); plot(snr,lcc,'k-'); grid on; ylim([1,10]);
hold off
xlabel('SNR (dB)'); ylabel('Capacity (bit/symbol)'); legend('QPSK','8-QAM','16-QAM','32-QAM','64-QAM');


