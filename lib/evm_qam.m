function evm = evm_qam(x)
x = DspAlg.Normalize(x,16);
% decision
bound = 2;
rx = [sign(real(x)), sign(real(x)-sign(real(x))*bound), ...
      sign(imag(x)), sign(imag(x)-sign(imag(x))*bound)];
dx = sign(rx+1);
% bit 2 sym
s = [dx(:,1)*8, dx(:,2)*4, dx(:,3)*2, dx(:,4)];
sx = sum(s,2) + 1;

for ii = 1:16
    idx{ii} = find(sx==ii);
    sig(ii) = mean(x(idx{ii}));
    n(ii) = mean(abs(x(idx{ii})-sig(ii)).^2);
end

evm = sum(n) / sum(abs(sig).^2);