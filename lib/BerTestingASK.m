function bn = BerTestingASK(signal, mn, BitReference, differential)
x = dsp.Normalize(signal, mn);
Rx = MakeDecision(x, mn);
if differential
    Rx = DiffDecodeBinary(Rx);
    Tx = DiffDecodeBinary(BitReference);
    bn = CompareRxTx(Rx, Tx);
else
    Tx1 = BitReference;
    Tx2 = ~BitReference;
    bn = min( CompareRxTx(Rx, Tx1),CompareRxTx(Rx, Tx2) );
end
return

function y = DiffDecodeBinary(x)
x_delay = [0; x(1:end-1)];
y = xor(x, x_delay);
return

function BN = CompareRxTx(Rx, Tx)
% set the window size to 1000
n = 1000;    BN = 0;
% double the reference
LongTx = [Tx; Tx];
% creat a window from the received signal
refpattern = Rx(1:n);
% searching...
for k = 1:length(Tx)
    datpattern = LongTx(k:k+n-1);
    bn_tmp(k) = nnz(datpattern - refpattern);
end
[mbn, idx] = min(bn_tmp);
% figure; plot(bn_tmp,'.-'); ylim([0 n]);
% legend(['Find minimun error ' num2str(mbn) 'at idx ' num2str(idx)])
% comparing...
if length(Rx) <= length(Tx)
    ref = LongTx(idx:idx+length(Rx)-1);
    dat = Rx;
    BN = nnz(dat - ref);
else
    leng2 = mod(length(Rx), length(Tx));
    leng1 = length(Rx) - leng2;
    ref = LongTx(idx:idx+length(Tx)-1);
    for k = 1:length(Tx):leng1
        dat = Rx(k:k+length(Tx)-1);
        BN = BN + nnz(dat - ref);
    end
    ref2 = ref(1:leng2);
    dat2 = Rx(leng1+1:end);
    BN = BN + nnz(dat2-ref2);
end
return

function b = MakeDecision(x, mn)
switch mn
    case 2
        rx = sign(real(x));
    case 1
        rx = sign(abs(x) - mean(abs(x)));
end
b = sign( rx + 1 );
return
