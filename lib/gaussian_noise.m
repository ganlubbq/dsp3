function w = gaussian_noise(nrow, ncol, p, pmod, dmod)
% Generate real or complex white Gaussian noise
% 
% w = gaussian_noise(nrow, ncol, p, pmod, dmod) 
%   p : noise power or variance

% default impedance
imp = 1;
w = zeros(nrow, ncol);

switch lower(pmod)
    case 'linear'
        noisePower = p;
    case 'dbw'
        noisePower = 10.^(.1 * p);
    case 'dbm'
        noisePower = 10.^(.1 * (p - 30));
    otherwise
        warning('noise::power mode'); keyboard;
end

switch lower(dmod)
    case 'complex'
        if nrow > ncol
            for ii = 1 : ncol
                rnd0 = randn(nrow, 1); rnd0n = rnd0 / sqrt(mean(rnd0.^2));
                rnd1 = randn(nrow, 1); rnd1n = rnd1 / sqrt(mean(rnd1.^2));
                w(:, ii) = (sqrt(imp * noisePower / 2)) * (rnd0n + 1i * rnd1n);
            end
        else
            for ii = 1 : nrow
                rnd0 = randn(1, ncol); rnd0n = rnd0 / sqrt(mean(rnd0.^2));
                rnd1 = randn(1, ncol); rnd1n = rnd1 / sqrt(mean(rnd1.^2));
                w(ii, :) = (sqrt(imp * noisePower / 2)) * (rnd0n + 1i * rnd1n);
            end
        end
    case 'real'
        if nrow > ncol
            for ii = 1 : ncol
                rnd0 = randn(nrow, 1); rnd0n = rnd0 / sqrt(mean(rnd0.^2));
                w(:, ii) = (sqrt(imp * noisePower)) * rnd0n;
            end
        else
            for ii = 1 : nrow
                rnd0 = randn(1, ncol); rnd0n = rnd0 / sqrt(mean(rnd0.^2));
                w(ii, :) = (sqrt(imp * noisePower)) * rnd0n;
            end
        end
    otherwise
        warning('noise::data mode'); keyboard;
end
return
