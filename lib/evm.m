function v = evm(x, mn)
% Get error vector magnitude for input signals

xn = normalizeQam(x, mn);
xd = slicerGrayQam(xn, mn);
Ps = mean(abs(xd).^2);
c = constellation(mn);

for ii = 1 : mn
    idx{ii} = find(xd == c(ii));
    n(ii) = mean(abs(xn(idx{ii}) - c(ii)).^2);
end

v = mean(n) / Ps;
% v = sqrt(v);

return
