function phase_noise = phase_noise(nsample, pnvar, p_ini)
% Generate random walk phase noise

tmp = randn(1, nsample);

% remove dc component
tmp = tmp - mean(tmp);

% normalize
tmp = tmp ./ calcrms(tmp);
phase_noise = p_ini + cumsum(tmp .* sqrt(pnvar));
return
