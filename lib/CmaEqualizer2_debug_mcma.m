function [Y h1 h2 df] = CmaEqualizer2_debug_mcma(xx, h1, h2, taps, mu, R, sps, err_id, stage, method)
%CMAADAPTIVEFILTER Polarization demultiplexing filter using CMA algorithm
%
%   Example
%   
%   See also PolarizationDemux,CmaEqualizer

%   Copyright2010 WANGDAWEI 16/3/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 10, method = 1;     end
if nargin < 9,  stage = 2;      end
if nargin < 8,  err_id = 1;     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(xx);             % Length of eXtended X = length(X) + taps -1
ntap = taps;                % Number of Taps
L = N-ntap+1;               % Length of outputs
Y1 = zeros( L ,1 );         % First output vector initialized to zero
Y2 = zeros( L ,1 );         % Second output vector initialized to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dontskip = (sps == 1);      % sps>1, updating only once per sps samples
modk = 0;                   % updating phase [0 ~ sps-1]
                            % Note that the updating phase determines the downsamping phase
                            % of the equalizer output in function 'PolarizationDemux'
df = zeros( L ,1 );
s4 = 0;
Y1(1) = 1+1j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ktime = 2:L
    nindex = ktime : ktime+ntap-1;      % time index for inputs
    X1 = xx( nindex, : ) .* h1;
    X2 = xx( nindex, : ) .* h2;
    Y1( ktime ) =  sum( sum( X1 ) );
    Y2( ktime ) =  sum( sum( X2 ) );
    if stage == 1                       % 1st stage CMA and update h1 only
        if (dontskip || mod(ktime-1,sps)==modk )
            incr1 = mu.*errorfuncma( Y1(ktime), R, err_id).*conj( xx(nindex,:) );
            h1 = h1 + incr1;
            h2(:,2) = ConjReverse( h1(:,1) );
            h2(:,1) = -1*ConjReverse( h1(:,2) );
        end
    else                                % 2nd stage CMA and update h1h2
        if (dontskip || mod(ktime-1,sps)==modk)
            incr1 = mu.*errorfuncma( Y1(ktime), R, err_id).*conj( xx(nindex,:) );
            incr2 = mu.*errorfuncma( Y2(ktime), R, err_id).*conj( xx(nindex,:) );
            h1 = h1 + incr1;
            h2 = h2 + incr2;
        end
        s4 = s4 + Y1(ktime).^4*conj(Y1(ktime-1).^4);
        df(ktime+1) =  -angle(s4)/4;
    end
    switch method
    end
end
fo = (1:L)*mean(df);
Y = [Y1.*exp(1j*fo(:)) Y2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Err = errorfuncma( X, R, IDX)
switch IDX
    case 1   %%% Classical CMA
        Err = -X.*(abs(X).^2-R(1));
    case 2   %%% Modified CMA
        Err = -complex(real(X).*(real(X).^2-R(1)),imag(X).*(imag(X).^2-R(1)));
    case 3   %%% Cascaded multi-modulus algorithm
        A1 = 0.5*(R(1)+R(2));
        A2 = 0.5*(R(3)-R(1));
        A3 = 0.5*(R(3)-R(2));
        e1 = abs(X)-A1; e2 = abs(e1)-A2; e3 = abs(e2)-A3;
        Err = -sign(X) * e3 * sign(e1) * sign(e2);
    case 4   %%% Radius directed
        Err = -X.*min(abs(abs(X).^2-R));
    case 5
        Err = complex(real(X).*(10 - 2*real(X).^2*9),imag(X).*(10 - 2*imag(X).^2*9));
    case 6
        Err = complex(real(X).*(4 - abs(real(X)^2*9-5)),imag(X).*(4 - abs(imag(X)^2*9-5)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ConjReverse(x)
y = conj(flipud(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%