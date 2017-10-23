function v = evm(x, mn)
% Get error vector magnitude for input signals

xn = normalization(x, mn);
xd = slicer_mqam(xn, mn);
ps = mean(abs(xd).^2);
cc = constellation(mn);

for ii = 1 : mn
    idx{ii} = find(xd == cc(ii));
    n(ii) = mean(abs(xn(idx{ii}) - cc(ii)).^2);
end

v = mean(n) / ps;
% v = sqrt(v);
return
