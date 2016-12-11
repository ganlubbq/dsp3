% Generate real or complex white Gaussian noise
% 
% Example: y = genWGN(nrow,ncol,p,pmod,dmod)
% 
% Input: 
%       len     - length of data
%       p       - noise power or variance
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function y = genWGN(nrow,ncol,p,pmod,dmod)

% default impedance
imp = 1;
y = zeros(nrow,ncol);

switch lower(pmod)
    case 'linear'
        noisePower = p;
    case 'dbw'
        noisePower = 10^(p/10);
    case 'dbm'
        noisePower = 10^((p-30)/10);
    otherwise
        error('incorrect noise power mode');
end
ss = nrow>ncol;
switch lower(dmod)
    case 'complex'
        if ss
            for ii = 1:ncol
                rnd0 = randn(nrow,1); rnd0n = rnd0/sqrt(mean(rnd0.^2));
                rnd1 = randn(nrow,1); rnd1n = rnd1/sqrt(mean(rnd1.^2));
                y(:,ii) = (sqrt(imp*noisePower/2))*(rnd0n+1j*rnd1n);
            end
        else
            for ii = 1:nrow
                rnd0 = randn(1,ncol); rnd0n = rnd0/sqrt(mean(rnd0.^2));
                rnd1 = randn(1,ncol); rnd1n = rnd1/sqrt(mean(rnd1.^2));
                y(ii,:) = (sqrt(imp*noisePower/2))*(rnd0n+1j*rnd1n);
            end
        end
    case 'real'
        if ss
            for ii = 1:ncol
                rnd0 = randn(nrow,1); rnd0n = rnd0/sqrt(mean(rnd0.^2));
                y(:,ii) = (sqrt(imp*noisePower))*rnd0n;
            end
        else
            for ii = 1:nrow
                rnd0 = randn(1,ncol); rnd0n = rnd0/sqrt(mean(rnd0.^2));
                y(ii,:) = (sqrt(imp*noisePower))*rnd0n;
            end
        end
    otherwise
        error('Incorrect date mode in genWGN!!!');
end

return


