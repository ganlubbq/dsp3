

function [ER,EC] = matlab_bert(signal,mn,reference,method)
%BERTESTING Count the number of bit error
%
%   Example
%   
%   See also: CoderDriver

%   Copyright2010 WANGDAWEI $16/3/2010$

cutoff = 8;

signal(isnan(signal)) = 1+1i;
signal(signal==0) = 1+1i;

% normalize the data
tmp_x = DspAlg.Normalize(signal,mn);

% loading reference
tmp_Tx1 = reference{1};
tmp_Tx2 = reference{2};

% cut-off
x   = tmp_x(1+cutoff:end-cutoff,:);
Tx1 = tmp_Tx1(1+cutoff:end-cutoff);
Tx2 = tmp_Tx2(1+cutoff:end-cutoff);

if nargin<4
    method = 'diff';
end

if strcmp(method,'binary') || strcmp(method,'gray')
    Tx_bits1 = matlab_demod(Tx1,mn,method);
    Tx_bits2 = matlab_demod(Tx2,mn,method);
    Tx_bits1 = Tx_bits1.';
    Tx_bits2 = Tx_bits2.';
    rk = 3;
elseif strcmp(method,'diff')
    Tx_bits1 = matlab_demod(Tx1,mn,'binary');
    Tx_bits2 = matlab_demod(Tx2,mn,'binary');
    Tx_bits1 = DspAlg.DifferentialDecode(Tx_bits1.', mn);
    Tx_bits2 = DspAlg.DifferentialDecode(Tx_bits2.', mn);
    rk = 0;
else
    error('encode/decode method not supported.');
end

for rotate = 0:rk
    % rotate and try
    x = x.*exp(-1j*rotate*pi/2);
    for pol = 1:2
        if strcmp(method,'binary') || strcmp(method,'gray')
            Rx_bits = matlab_demod(x(:,pol),mn,method);
            Rx_bits = reshape(Rx_bits,log2(mn),[]);
            Rx_bits = Rx_bits.';
        elseif strcmp(method,'diff')
            Rx_bits = matlab_demod(x(:,pol),mn,'binary');
            Rx_bits = reshape(Rx_bits,log2(mn),[]);
            Rx_bits = DspAlg.DifferentialDecode(Rx_bits.', mn);
        end
        
        % we dont know which polar is which, test.
        n1 = CompareRxTx(Rx_bits,Tx_bits1);
        n2 = CompareRxTx(Rx_bits,Tx_bits2);
        
        EC(rotate+1,pol) = (n1 <= n2)*n1 + (n1 > n2)*n2;
    end
end

[~,idx_rot] = min(EC,[],1);
idx_rot = idx_rot-1;

% ec should be doubled if signal is differential decoded
EC = sum(min(EC,[],1));
ER = EC/(2*length(x)*log2(mn));


function n = CompareRxTx(Rx, Tx)
% a small shift length for comparison
Reflen = 1000;
% search
for k = 1 : 10
    Rxy(k) = nnz(xor(Rx(1:1+Reflen,:),Tx(k:k+Reflen,:)));
    Ryx(k) = nnz(xor(Tx(1:1+Reflen,:),Rx(k:k+Reflen,:)));
end
[val_ahead idx_ahead] = min(Rxy);
[val_delay idx_delay] = min(Ryx);
if val_ahead < val_delay
    inx = idx_ahead;
    Ntest = min(length(Tx)-inx +1, length(Rx));
    n = nnz(Tx(inx:inx+Ntest-1,:)-Rx(1:Ntest,:));
else
    inx = idx_delay;
    Ntest = min(length(Rx)-inx +1, length(Tx));
    n = nnz(Rx(inx:inx+Ntest-1,:)-Tx(1:Ntest,:));
end

