function [y phi] = cpe_kalman(x,mn,mu,appML,bs,iter)
%FEEDFACKCPE Decision-directed feefback carrier phase recovery routine
%
% Example
% 
% See also FeedforwardCPE
%
% [1] X. Yang, X. Cui, M. Lu, and Z. Fen, ¡°Carrier recovery using FFT and
% Kalman filter,¡± in International Symposium on Image and Signal Processing
% and Analysis, 2003. ISPA 2003., 2003, vol. 2, pp. 1094¨C1096.

% Copyright2010 WANGDAWEI 16/3/2010

if nargin<6
    iter = 1;
end
if nargin<5
    bs = 16;
end
if nargin<4
    appML = 0;
end
x = DspAlg.Normalize(x,mn);
x = [x(1,:);x];

y = zeros(size(x));
phi = zeros(size(x));
err = zeros(size(x));

N = length(x);

for kk = 2:N
    err(kk,:) = mod(angle(x(kk,:))-phi(kk-1,:),pi/4);
    phi(kk,:) = phi(kk-1,:) + mu*err(kk,:);
    y(kk,:) = x(kk,:).* exp(-1j*phi(kk,:));
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
