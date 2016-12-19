function [ER EC] = BerTestingQAM(signal, mn, reference, diffFlag)

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
% set the window size to 1000
n = 1000;
% double the reference
LongTx = [Tx; Tx];
% creat a window from the received signal
refpattern = Rx(1:n,:);
% searching...
for k = 1:length(Tx)
    datpattern = LongTx(k:k+n-1,:);
    bn_tmp(k) = nnz(datpattern - refpattern);
end
[~,idx] = min(bn_tmp);
% comparing...
ref = LongTx(idx:length(Tx),:);
dat = Rx(1:length(Tx)-idx+1,:);
BN = nnz(dat - ref);
return
