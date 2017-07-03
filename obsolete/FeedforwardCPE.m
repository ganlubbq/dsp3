function [y,phn] = FeedforwardCPE(x,mn,bs,method,appML,iter)
%FEEDFORWARDCPE Feedforward carrier phase estimation routine.
%   bs: block size
%   method: block,slide.
%   iter: iteration of ML stage
%
%   VVPE with qpsk partition including both block-by-block and slide window
%   method for 16-QAM signal.
%
%   Example
%   
%   See also FeedbackCPE

%   copyright2010 wangdawei 16/3/2010

if nargin<6
    iter = 1;
end
if nargin<5
    appML = 0;
end
if nargin<4
    method = 'block';
end

x(isnan(x)) = 1+1i;

% add normalize here
x = DspAlg.Normalize(x,mn);

if strcmp(method,'bps') || mn == 64
    [y phn] = bps_cpe(x,bs,90,mn);
elseif mn==2
    % raise the phase to 2nd power not the amplitude
    x2 = abs(x).*exp( 1j*angle( x.^ 2 ) );
    % sliding average
    for pol = 1:size(x,2)
        x2(:,pol) = smooth(x2(:,pol),bs);
    end
    % unwrapping
    phn = unwrap(angle(x2))/2;
    % rotating
    y = x.* exp(-1j*phn);
elseif mn==4
    % raise the phase to 4th power not the amplitude
    x4 = abs(x).*exp( 1j*angle( x.^ 4 ) );
    % sliding average
    for pol = 1:size(x,2)
        x4(:,pol) = smooth(x4(:,pol),bs);
    end
    % unwrapping
    phn = unwrap(angle(x4))/4;
    % rotating
    y = x.* exp(-1j* (phn + pi/4) );
elseif mn==16
    if strcmp(method,'block')
        [y phn] = block_cpe(x,bs);
    end
    if strcmp(method,'slide')
        [y phn] = slide_cpe(x,bs);
    end
end

if appML
    b = zeros(size(y));
    h = zeros(size(y));
    for ii = 1:iter
        for pol = 1:size(y,2)
            b(:,pol) = DspAlg.slicer(y(:,pol),mn);
            h(:,pol) = smooth(y(:,pol).*conj(b(:,pol)),bs);
            y(:,pol) = y(:,pol).* exp(-1j*angle(h(:,pol)));
        end
    end
end


function [y phn] = bps_cpe(x,bs,M,mn)
% How to solve the cycle-slip problem
phi = (0:M-1) / M * pi/2;
y = zeros(size(x));

% distance
D = zeros(size(x));

% phase noise
phn = zeros(size(x));

for pol = 1:size(x,2)
    xp = repmat(x(:,pol),1,M);
    xd = zeros(size(xp));
    for p = 1:M
        xp(:,p) = xp(:,p).* exp(1j*phi(p));
        xd(:,p) = DspAlg.slicer(xp(:,p),mn);
        D(:,p) = smooth(abs(xp(:,p)-xd(:,p)).^2,bs);
    end
    [~,idx] = min(D,[],2);
    rot_phi = phi(idx);
    rot_phi_un = unwrap(4*rot_phi)/4;
    y(:,pol) = x(:,pol).* exp(1j*rot_phi_un(:));
    phn(:,pol) = rot_phi_un(:);
end


function [y phn] = block_cpe(x, bs)
% 
a = mod(size(x,1),bs);
if a
    x = [x;zeros(bs-a,size(x,2))];
end

% two reference circles
R1 = (sqrt(2)+sqrt(10))/2;
R2 = (3*sqrt(2)+sqrt(10))/2;

for pol = 1:size(x,2)
    
    x_bk = reshape(x(:,pol),bs,[]);
    
    idx = abs(x_bk)<R1 | abs(x_bk)>R2;
    d1 = x_bk.* idx;
    s1 = ones(bs,1) * sum(sign(d1.^4));
    
    d2 = x_bk.* ~idx;
    r1 = sign((d2*exp(1j*(pi/4-atan(1/3)))).^4);
    r2 = sign((d2*exp(1j*(atan(1/3)-pi/4))).^4);
    idx = r1~=0 & (abs(r1 - sign(s1))<= abs(r2 - sign(s1)));
    d3 = r1.* idx;
    d4 = r2.*~idx;
    s2 = sum(d3) + sum(d4);
    
    sume = s1(1,:) + s2;
    
    PN  = unwrap(angle(sume)+1*pi)./4;
    phi = repmat(PN, bs, 1);
    x_bk_com = x_bk.*exp(-1j*phi);
    y(:,pol) = x_bk_com(1:end-(a>0)*(bs-a));
    phn(:,pol) = phi(1:end-(a>0)*(bs-a));
end


function [y, phn] = slide_cpe(x, bs)
bs = bs + ~mod(bs,2);
y = zeros(size(x));
[nrow, ncol] = size(y);
halfN = floor(bs/2);
xx = [x(end-halfN+1:end,:); x; x(1:halfN,:)];
for ii = 1 : ncol
    x = xx(:,ii);
    sum_1 = 0;
    sum_2 = 0;
    sum_sum = zeros(size(x));
    % partition and orientate
    for k = halfN+1:nrow+halfN
        for kk = k-halfN:k+halfN
            if abs(x(kk))<(sqrt(2)+sqrt(10))/2 || abs(x(kk))>(3*sqrt(2)+sqrt(10))/2
                sum_1 = sum_1 + sign(x(kk).^4);
            else
                rdata1 = sign((x(kk)*exp(1j*(pi/4-atan(1/3)))).^4);
                rdata2 = sign((x(kk)*exp(1j*(atan(1/3)-pi/4))).^4);
                if abs(rdata1 - sign(sum_1)) <= abs(rdata2 - sign(sum_1))
                    sum_2 = sum_2 + 2*rdata1;
                else
                    sum_2 = sum_2 + 2*rdata2;
                end
            end
        end
        sum_sum(k) = (sum_1 + sum_2)/bs;
        sum_1 = 0; sum_2 = 0;
    end
    phn = unwrap( angle(sum_sum(halfN+1:end-halfN))+pi )/4;
    y(:,ii) = x(halfN+1:end-halfN).*exp(-1j*phn);
end
