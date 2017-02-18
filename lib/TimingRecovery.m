% DESCRIPTION
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright default

function [ y, tphase ] = TimingRecovery( signal, mn, sps, mu, method )

switch method
    case 'simple'
        y = signal(1:sps:end, :);
        tphase = [0;0];
    case {'linear','parabolic'}
        x = dsp.Normalize(signal, mn);
        [y,tphase] = GardnerPLL(x, sps, mu, method);
    otherwise
        error('unsupported timing recovery method !!!');
end

return

function [symRtm phEst]= GardnerPLL(x, sps, g, method)
T  = 1;                      % Normalized symbol period
Ts = 0.5;                    % Half a symbol period
N = sps;
Ncol = size(x, 2);
m = mod(length(x), N);
if m,  x = [x; zeros(N-m, Ncol)];  end
numSym = length(x) / N;
numBuf = (numSym + 3) * N;   % +3 symbols overlap
buffer = zeros(numBuf,1);
phEst  = zeros(numSym,Ncol);
symRtm = zeros(numSym,Ncol);
for jj = 1 : Ncol
tau = zeros(numSym+2, 1);
buffer( 1+N : numBuf-N-N ) = x(:,jj);
if strcmp(method, 'linear')
    for k = 1 : numSym
        % Sampling instants for tau calculation
        fltIdx1 = N*( (k-1)*T  + tau(k));
        fltIdx2 = N*( k*T      + tau(k+1));
        fltIdx3 = N*( k*T - Ts + tau(k));
        % Select the second symbol as the current symbol so as to have...
        % a buffer of symbols on both sides [-T..T..+2T]
        idx1_dn = N+1 + floor(fltIdx1);
        idx2_dn = N+1 + floor(fltIdx2);
        idx3_dn = N+1 + floor(fltIdx3);
        % interpolate between 2 available points
        a1 = linearInterp(buffer, idx1_dn, fltIdx1);
        a2 = linearInterp(buffer, idx2_dn, fltIdx2);
        a3 = linearInterp(buffer, idx3_dn, fltIdx3);
        % Timing error calculation
        b1 = real(a1 - a2) * real(a3);
        b2 = imag(a1 - a2) * imag(a3);
        eK = b1 + b2;
        % Apply current estimate to next. Actrually we need k and k+1..
        % errors to calculate current error and output current symbol...
        % so k+2 error should be updated before we continue next iteration
        tau(k+2) = tau(k+1) + g * eK;
        if abs(tau(k+2)) > Ts
            tau(k+2) = rem(tau(k+2) + Ts, T); % [-T, T], due to fmod */
            if tau(k+2) < 0
                tau(k+2) = tau(k+2)+ T;   % [0, T] */
            end
            tau(k+2) = tau(k+2)-Ts;       % [-T/2, T/2] */
        end
    end
    for k = 1 : numSym
        % always delay so tau must > zero
        if tau(k) < 0
            tauInt = tau(k) + T;
        else
            tauInt = tau(k);
        end
        % output symbol phase estimate
        phEst(k,jj) = tauInt * N;
        dn = k*N + 1 + floor(phEst(k,jj));
        symRtm(k,jj) = linearInterp(buffer, dn, phEst(k,jj));
    end
elseif strcmp(method, 'parabolic')
    for k = 1 : numSym
        % Sampling instants for tau calculation
        fltIdx1 = N*( (k-1)*T  + tau(k));
        fltIdx2 = N*( k*T      + tau(k+1));
        fltIdx3 = N*( k*T - Ts + tau(k));
        % Select the second symbol as the current symbol so as to have...
        % a buffer of symbols on both sides [-T..T..+2T]
        idx1_dn = N+1 + floor(fltIdx1);
        idx2_dn = N+1 + floor(fltIdx2);
        idx3_dn = N+1 + floor(fltIdx3);
        % interpolate between 2 available points
        a1 = ParabolicInterp(buffer, idx1_dn, fltIdx1);
        a2 = ParabolicInterp(buffer, idx2_dn, fltIdx2);
        a3 = ParabolicInterp(buffer, idx3_dn, fltIdx3);
        % Timing error calculation
        b1 = real(a1 - a2) * real(a3);
        b2 = imag(a1 - a2) * imag(a3);
        eK = b1 + b2;
        % Apply current estimate to next. Actrually we need k and k+1..
        % errors to calculate current error and output current symbol...
        % so k+2 error should be updated before we continue next iteration
        tau(k+2) = tau(k+1) + g * eK;
        if abs(tau(k+2)) > Ts
            tau(k+2) = rem(tau(k+2) + Ts, T); % [-T, T], due to fmod */
            if tau(k+2) < 0
                tau(k+2) = tau(k+2)+ T;   % [0, T] */
            end
            tau(k+2) = tau(k+2)-Ts;       % [-T/2, T/2] */
        end
    end
    for k = 1 : numSym
        % always delay so tau must > zero
        if tau(k) < 0
            tauInt = tau(k) + T;
        else
            tauInt = tau(k);
        end
        % output symbol phase estimate
        phEst(k,jj) = tauInt * N;
        dn = k*N + 1 + floor(phEst(k,jj));
        symRtm(k,jj) = ParabolicInterp(buffer, dn, phEst(k,jj));
    end
end
return


function y = linearInterp(buffer, dn, fltIdx)
a = real(buffer(dn)) + (fltIdx - floor(fltIdx)) * (real(buffer(dn+1)) - real(buffer(dn)));
b = imag(buffer(dn)) + (fltIdx - floor(fltIdx)) * (imag(buffer(dn+1)) - imag(buffer(dn)));
y = a + 1j*b;
return


function y = ParabolicInterp(buffer, dn, fltIdx)
u = fltIdx-floor(fltIdx);
a = 0.5;
C1 = a*u*u - a*u;
C0 = -1*a*u*u + (a-1)*u + 1;
C_1 = -1*a*u*u + (a+1)*u;
C_2 = a*u*u - a*u;
a = C1 *real(buffer(dn-1)) + C0 * real(buffer(dn)) + C_1*real(buffer(dn+1)) + C_2*real(buffer(dn+2));
b = C1 *imag(buffer(dn-1)) + C0 * imag(buffer(dn)) + C_1*imag(buffer(dn+1)) + C_2*imag(buffer(dn+2));
y = a + 1j*b;
return

