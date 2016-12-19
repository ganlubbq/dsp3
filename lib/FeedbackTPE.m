function [y,phEst,err] = FeedbackTPE(signal,mn,sps,mu,estMeth,intMeth,decFlag,norFlag)
%TIMINGRECOVERY Summary of this function goes here
%   mn: modulation order
%   mu: step size
%   estMeth: estimator method
%   intMeth: interpolation method
%
%   Example
%   
%   See also FeedforwardTPE

%   Copyright2010 WANGDAWEI 16/3/2010

if nargin<8
    norFlag = 1;
end
if nargin<7
    decFlag = 1;
end
if nargin<6
    intMeth = 'linear';
end
if nargin<5
    estMeth = 'gardner';
end

% normalize the input or not
if norFlag
    x = DspAlg.Normalize( signal, mn)/(sqrt(mn)-1);
else
    x = signal;
end

% get the size of input
mm = size(x,1);
kk = size(x,2);

% chose the estimator
switch lower(estMeth)
    case 'none'
        y       = x(1:sps:end,:);
        phEst   = zeros(1,kk);
        err     = zeros(1,kk);
    case 'gardner'
        % make sure that the length of x can be diveded by sps
        temp = mod(mm, sps);
        if temp
            x = [x;zeros(sps-temp,kk)];
        end
        for ii = 1:kk
            [y(:,ii),phEst(:,ii),err(:,ii)] = GardnerPLL(x(:,ii),sps,mu,intMeth);
        end
end
if decFlag
    y = y(1:2:end,:);
end

function [symRtm,phEst,ee] = GardnerPLL(x,sps,g,method)
%GARDNERPLL Gardner timing error detector
%
% The modified version is refered to:
%
% [1] W. Gappmair, S. Cioni, G. E. Corazza, and O. Koudelka, ��Symbol-Timing
% Recovery with Modified Gardner Detectors,�� in International Symposium on
% Wireless Communication Systems, 2005, pp. 831�C834.
%
% When mu = 1, the ted reduce to the normal gardner ted.
%

T   = 1;                        % Normalized symbol period
Ts  = 0.5;                      % Half a symbol period
N   = sps;

numSym  = length(x) / N;
numBuf  = ( numSym+3 ) * N;     % +3 symbols overlap
buffer  = zeros(numBuf,1);
phEst   = zeros(numSym,1);
symRtm  = zeros(numSym*N,1);
ee      = zeros(numSym,1);
tau     = zeros(numSym+2,1);

buffer(N+1:numBuf-N-N) = x;

for k = 1 : numSym
    
    % ================= 1st part =================== %
    % Sampling instants for tau calculation
    fltIdx1 = N * tau(k);
    fltIdx3 = N * tau(k+1);
    fltIdx2 = N * tau(k);
    
    % Select the second symbol as the current symbol so as to have...
    % a buffer of symbols on both sides [-T..T..+2T]
    idx1_dn = k*N + 1 + floor(fltIdx1);
    idx3_dn = (k+1)*N + 1 + floor(fltIdx3);
    idx2_dn = (k+1/2)*N + 1 + floor(fltIdx2);
    
    % interpolate between 2 available points [interpolater output]
    a1 = Interp(buffer, idx1_dn, fltIdx1, method); % 1st sample of 1st symbol
    a3 = Interp(buffer, idx3_dn, fltIdx3, method); % 1st sample of 2nd symbol
    a2 = Interp(buffer, idx2_dn, fltIdx2, method); % 2nd sample of 1st symbol
    
    % the following output has slight difference with symRtm when tau<0
    symRtm(2*k-1) = a1;
    symRtm(2*k) = a2;
    
    % ================= 2nd part =================== %
    % Timing error detector
    b1 = real(a3 - a1) * real(a2);
    b2 = imag(a3 - a1) * imag(a2);
    eK = b1 + b2;
    
    % Modified timing error detector
    %mu = 1;
    %px3 = abs(a3).^ mu.* exp(1i*angle(a3));
    %px1 = abs(a1).^ mu.* exp(1i*angle(a1));
    %eK = real((px3-px1).* conj(a2));

    % ================= 3rd part =================== %
    % Apply current estimate to next. Actrually we need k and k+1..
    % errors to calculate current error and output current symbol...
    % so k+2 error should be updated before we continue next iteration
    tau(k+2) = tau(k+1) - g * eK;
    ee(k) = eK;
    
    % wrap timing phase
    if abs(tau(k+2)) > Ts
        tau(k+2) = rem(tau(k+2) + Ts, T);   % [-T, T], due to rem
        if tau(k+2) < 0
            tau(k+2) = tau(k+2)+ T;         % [0, T]
        end
        tau(k+2) = tau(k+2)-Ts;             % [-T/2, T/2]
    end
    
    phEst(k) = tau(k) * N;
end


function y = Interp(buffer,idx,epsilon,method)
%INTERP Interpolate signal samples
%

mu = epsilon - floor(epsilon);

if strcmpi(method,'linear')
    y = buffer(idx) + mu * (buffer(idx+1) - buffer(idx));
elseif strcmpi(method,'parabolic')
    a = 0.5;
    C_1 = a*mu*mu - a*mu;
    C0 = -a*mu*mu + (a-1)*mu + 1;
    C1 = -a*mu*mu + (a+1)*mu;
    C2 = a*mu*mu - a*mu;
    y = C_1 *buffer(idx-1) + C0 * buffer(idx) + C1*buffer(idx+1) + C2*buffer(idx+2);
end

