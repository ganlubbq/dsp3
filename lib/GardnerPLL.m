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
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function [symRtm, phEst]= GardnerPLL(symIn, sps, gain, int_method)

% Normalized symbol period
T  = 1;
% Half a symbol period
Ts = 0.5;
N = sps;

% process along data columns
Ncol = size(symIn, 2);

% padding zeros
m = mod(length(symIn), N);
if m,  symIn = [symIn; zeros(N-m, Ncol)];  end

numSym = length(symIn) / N;

numBuf = (numSym + 3) * N;   % +3 symbols overlap

% initializing
buffer = zeros(numBuf,1);
phEst  = zeros(numSym,Ncol);
symRtm = zeros(numSym,Ncol);

% process along data columns
for jj = 1 : Ncol
    
    % flush everything
    tau = zeros(numSym+2, 1);
    buffer( 1+N : numBuf-N-N ) = symIn(:,jj);
    
    if strcmp(int_method, 'linear')
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
            tau(k+2) = tau(k+1) + gain * eK;
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
    elseif strcmp(int_method, 'parabolic')
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
            tau(k+2) = tau(k+1) + gain * eK;
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
end
return

% linear interpolation using 2 samples
function y = linearInterp(buffer, dn, fltIdx)
a = real(buffer(dn)) + (fltIdx - floor(fltIdx)) * (real(buffer(dn+1)) - real(buffer(dn)));
b = imag(buffer(dn)) + (fltIdx - floor(fltIdx)) * (imag(buffer(dn+1)) - imag(buffer(dn)));
y = a + 1j*b;
return

% parabolic interpolation using 4 samples
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


