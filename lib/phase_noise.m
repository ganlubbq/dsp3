function phase_noise = phase_noise(nsample, pvar, pini)
% Generate random walk discrete phase noise
tmp = randn(1, nsample);
% remove dc component followed by normalization
tmp = tmp - mean(tmp);
tmp = tmp ./ calcrms(tmp);
phase_noise = pini + cumsum(tmp .* sqrt(pvar));
return
