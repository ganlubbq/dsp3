function [Y h1 h2 err] = CmaEqualizer2_debug(xx, h1, h2, taps, mu, R, sps, err_id, stage, method)
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
halftaps = floor(taps/2);   % Number of half taps
L = N-ntap+1;               % Length of outputs
Y1 = zeros( L ,1 );         % First output vector initialized to zero
Y2 = zeros( L ,1 );         % Second output vector initialized to zero
err = zeros( 1 ,L );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dontskip = (sps == 1);      % sps>1, updating only once per sps samples
modk = 0;                   % updating phase [0 ~ sps-1]
                            % Note that the updating phase determines the downsamping phase
                            % of the equalizer output in function 'PolarizationDemux'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ktime = 1:L
    nindex = ktime : ktime+ntap-1;      % time index for inputs
    X1 = xx( nindex, : ) .* h1;
    X2 = xx( nindex, : ) .* h2;
    Y1( ktime ) =  sum( sum( X1 ) );    % Calculating outputs
    Y2( ktime ) =  sum( sum( X2 ) );
    err(ktime) = errorcalc( Y1(ktime), R, err_id );
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
    end
    switch method
    end
end
Y = [Y1 Y2];
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
function err = errorcalc(x,R,idx)
switch idx
    case 1   %%% Classical CMA
        err = -abs(x).^2-R(1);
    case 2   %%% Modified CMA
        err = -imag(x).^2-R(1);
    case 3   %%% Cascaded multi-modulus algorithm
        A1 = 0.5*(R(1)+R(2));
        A2 = 0.5*(R(3)-R(1));
        A3 = 0.5*(R(3)-R(2));
        e1 = abs(x)-A1; e2 = abs(e1)-A2; e3 = abs(e2)-A3;
        err = -e3 * sign(e1) * sign(e2) / abs(x);
    case 4   %%% Radius directed
        err = -min(abs(abs(x).^2-R));
    case 5
        err = 2*real(x).^2*9-10;
    case 6
        err = abs(real(x)^2*9-5)-4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ConjReverse(x)
y = conj(flipud(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%