% TEST SCRIPT FOR CALCULATING CONSTRAINED CHANNEL CAPACITY OF A DISCRETE
% MEMORYLESS GAUSSIAN CHANNEL
%
% Ungerboeck G. Channel coding with multilevel/phase signals[J]. 
% IEEE Transactions on Information Theory,
% 1982, 28(1):55-67.
%
%   ConstrainedCapacity([2 3])
%   ConstrainedCapacity([2 3 4])
%   ConstrainedCapacity([2 3 4 5])
%   ConstrainedCapacity([2 3 4 5 6])
%
% m-QAM
function ConstrainedCapacity(bitpersym)

if nargin < 1
    bitpersym = 2;
end

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',199));

symlen = 2^13;

for qq = 1 : length(bitpersym)
    mn = 2 .^ bitpersym(qq);
    
    % get the constellation
    a = constellation(mn);
    
    % input power
    sp = sum(abs(a).^2)/mn;
    
    % set snr range
    snr = 0 : 1 : 36;
    
    for pp = 1:length(snr)
        % get the power of noise
        th2 = 10*log10(sp) - snr(pp);
        
        % _TWO DIMENSION COMPLEX GAUSSIAN NOISE_
        z = wgn(1, symlen, th2, 'dbw', 'complex');
        
        for kk = 1:mn
            for ii = 1:mn
                tmp(ii,:) = exp(-(abs(a(kk) + z - a(ii)).^2 - abs(z).^2) / 10^(th2/10));
            end
            ss(kk,:) = log2(sum(tmp));
        end
        % the second term
        c2 = sum(mean(ss,2));
        
        chc(pp,qq) = log2(mn) - c2/mn;
    end
end

% Limit for complex modulation
lcc = log2(1 + idbw(snr));

figure(67);
hold on; 
plot(snr,lcc,'k-','LineWidth',2); 
plot(snr,chc,'LineWidth',2); 
grid on
ylim([1,10]);
xlabel('SNR (dB)'); 
ylabel('Capacity (bit/symbol)'); 
legend('Gaussian','QPSK','8-QAM','16-QAM','32-QAM','64-QAM');

return
