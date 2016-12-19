function evm = evm_qpsk(x)
% decision
dx = [sign(real(x)) sign(imag(x))];
dx = sign(dx+1);
% bit 2 sym
sx = dx(:,1)*2+dx(:,2);

idx0 = find(sx==0);
idx1 = find(sx==1);
idx2 = find(sx==2);
idx3 = find(sx==3);

sig0 = mean(x(idx0));
sig1 = mean(x(idx1));
sig2 = mean(x(idx2));
sig3 = mean(x(idx3));

n0 = mean(abs(x(idx0)-sig0).^2);
n1 = mean(abs(x(idx1)-sig1).^2);
n2 = mean(abs(x(idx2)-sig2).^2);
n3 = mean(abs(x(idx3)-sig3).^2);

evm = (n0+n1+n2+n3)/(abs(sig0)^2+abs(sig1)^2+abs(sig2)^2+abs(sig3)^2);