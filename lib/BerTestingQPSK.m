function [ER EC] = BerTestingQPSK(signal,mn,reference,diffFlag)

x = DspAlg.Normalize(signal, mn);
x(isnan(x)) = 0;

Rx1 = DspAlg.MakeDecision(x(:,1), mn);
Rx2 = DspAlg.MakeDecision(x(:,2), mn);

Tx1 = reference;

if diffFlag
    Rx1 = DspAlg.DifferentialDecode(Rx1, mn);
    Tx1 = DspAlg.DifferentialDecode(Tx1, mn);
    Rx2 = DspAlg.DifferentialDecode(Rx2, mn);
else
    Rx1 = DspAlg.sym2bit(DspAlg.bit2sym(Rx1,mn),mn);
    Tx1 = DspAlg.sym2bit(DspAlg.bit2sym(Tx1,mn),mn);
    Rx2 = DspAlg.sym2bit(DspAlg.bit2sym(Rx2,mn),mn);
end

EC = CompareRxTx(Rx1,Tx1) + CompareRxTx(Rx2,Tx1);

ER = EC/(2*length(x)*log2(mn));

return

function BN = CompareRxTx(Rx, Tx)
% periodic length
L = (2^7 - 1)*8;
% set the window size to 1000
n = 500;    BN = 0;
% double the reference
LongTx = [Tx(1:L,:); Tx(1:L,:)];
% creat a window from the received signal
refpattern = Rx(1:n,:);
% searching...
for k = 1 : L
    datpattern = LongTx(k:k+n-1,:);
    bn_tmp(k) = nnz(datpattern - refpattern);
end
[~, idx] = min(bn_tmp);
% comparing...
if length(Rx) <= L
    ref = LongTx(idx:idx+length(Rx)-1,:);
    dat = Rx;
    BN = nnz(dat - ref);
else
    ref = LongTx(idx:idx+L-1,:);
    leng2 = mod(length(Rx), L);
    leng1 = length(Rx) - leng2;
    for k = 1:L:leng1
        dat = Rx(k:k+L-1,:);
        BN = BN + nnz(dat - ref);
    end
    ref2 = ref(1:leng2,:);
    dat2 = Rx(leng1+1:end,:);
    BN = BN + nnz(dat2 - ref2);
end
return
