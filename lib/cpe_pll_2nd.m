% Decision-directed feefback carrier phase recovery routine implementing
% the phase lock loop
%
% Note:
% Example:
% 
% See also: 
%
% Referene: I. Fatadin, D. Ives, and S. J. Savory, "Compensation of Frequency
% Offset for Differentially Encoded 16- and 64-QAM in the Presence of Laser
% Phase Noise," IEEE Photonics Technology Letters, vol. 22, no. 3, pp.
% 176-178, Feb. 2010.
%
% Copyright2010 WANGDAWEI 16/3/2010

function [y phi] = cpe_pll_2nd(x,mn,gamma,rho,appML,bs,iter)

if nargin<7
    iter = 1;
end
if nargin<6
    bs = 16;
end
if nargin<5
    appML = 0;
end
x = DspAlg.Normalize(x,mn);
x = [x(1,:);x];
y = zeros(size(x));
zeta = zeros(size(x));
phi = zeros(size(x));
err = zeros(size(x));

N = length(x);

for kk = 2:N
    zz = x(kk,:) .* exp(-1j*phi(kk,:));
    y(kk,:) = zz;
    aa = slicer(zz,mn);
    err(kk,:) = imag( conj(aa).* zz );
    zeta(kk,:) = zeta(kk-1,:) + gamma*(1+rho)*err(kk,:) - gamma*err(kk-1,:);
    phi(kk+1,:) = phi(kk,:) + zeta(kk,:);
end

y = y(2:end,:);

if appML
    b = zeros(size(y));
    h = zeros(size(y));
    for ii = 1:iter
        for pol = 1:size(y,2)
            b(:,pol) = slicer(y(:,pol),mn);
            h(:,pol) = smooth(y(:,pol).*conj(b(:,pol)),bs);
            y(:,pol) = y(:,pol).* exp(-1j*angle(h(:,pol)));
        end
    end
end

return

