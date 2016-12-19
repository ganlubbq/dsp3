function [ER EC] = BerTesting(signal,mn,reference,method)
%BERTESTING Count the number of bit error
%   Test both polarizations by default, return the bit error number and the
%   bit error ratio. Actrually, the signal should be binary encoded.
%   Otherwise, if the signal is gray encoded, Tx should not be re-mapped
%   and the differential decode should have gray decode first.
%
%   Example
%   
%   See also: CoderDriver

%   copyright2010 WANGDAWEI $16/3/2010$

x = DspAlg.Normalize(signal,mn);
x(isnan(x)) = 0;

% loading reference
Tx1 = reference{1};  Tx1 = Tx1.';
Tx2 = reference{2};  Tx2 = Tx2.';

if nargin<4
    method = 'diff';
end

% mapping reference
switch method
    case 'diff'
        Tx1 = DspAlg.DifferentialDecode(Tx1, mn);
        Tx2 = DspAlg.DifferentialDecode(Tx2, mn);
        rk = 0;
    case 'gray'
        Tx1 = DspAlg.sym2bit(DspAlg.bit2sym(Tx1,mn),mn);
        Tx2 = DspAlg.sym2bit(DspAlg.bit2sym(Tx2,mn),mn);
        rk = 3;
    otherwise
        rk = 3;
end


for rotate = 0:rk
    
    % rotate and try
    x = x.*exp(-1j*rotate*pi/2);
    
    for pol = 1:2
        
        % hard decision
        Rx = DspAlg.MakeDecision(x(:,pol), mn);
        
        % mapping
        switch method
            case 'diff'
                Rx = DspAlg.DifferentialDecode(Rx, mn);
            case 'gray'
                Rx = DspAlg.sym2bit(DspAlg.bit2sym(Rx,mn),mn);
        end
        
        % we dont know which polar is which, test.
        n1 = CompareRxTx(Rx,Tx1);
        n2 = CompareRxTx(Rx,Tx2);
        
        EC(rotate+1,pol) = (n1 <= n2)*n1 + (n1 > n2)*n2;
        
    end
    
end

EC = sum(min(EC,[],1));
ER = EC/(2*length(x)*log2(mn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = CompareRxTx(Rx, Tx)
% a small shift length for comparison
Reflen = 1000;
% search
for k = 1 : 10
    Rxy(k) = nnz(~xor(Rx(1:1+Reflen,:),Tx(k:k+Reflen,:)));
    Ryx(k) = nnz(~xor(Tx(1:1+Reflen,:),Rx(k:k+Reflen,:)));
end
[val_ahead idx_ahead] = max(Rxy);
[val_delay idx_delay] = max(Ryx);
if val_ahead > val_delay
    inx = idx_ahead;
    Ntest = min(length(Tx)-inx +1, length(Rx));
    n = nnz(Tx(inx:inx+Ntest-1,:)-Rx(1:Ntest,:));
else
    inx = idx_delay;
    Ntest = min(length(Rx)-inx +1, length(Tx));
    n = nnz(Rx(inx:inx+Ntest-1,:)-Tx(1:Ntest,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF