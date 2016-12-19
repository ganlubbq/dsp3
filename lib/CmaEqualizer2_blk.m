function [Y h1 h2 mse deth] = CmaEqualizer2_blk(xx,h1,h2,taps, mu, R, sps, err_id, stage, method)
%CMAADAPTIVEFILTER Polarization demultiplexing filter using CMA algorithm
%   with block-by-block configuration
%
%   Example
%   
%   See also PolarizationDemux,CmaEqualizer2,CmaEqualizer

%   Copyright2010 WANGDAWEI 16/3/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 10, method = 1;     end
if nargin < 9,  stage = 2;      end
if nargin < 8,  err_id = 1;     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(xx,1);             % Length of eXtended X = length(X) + taps -1
ntap = taps;                % Number of Taps
x1 = reshape(xx(:,1),ntap,[]);
x2 = reshape(xx(:,2),ntap,[]);
nBlk = size(x1,2);
Y1 = zeros( N ,1 );         % First output vector initialized to zero
Y2 = zeros( N ,1 );         % Second output vector initialized to zero
mse = zeros( 1 ,nBlk );
deth = zeros( 1 ,nBlk );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dontskip = (sps == 1);      % sps>1, updating only once per sps samples
modk = 0;                   % updating phase [0 ~ sps-1]
                            % Note that the updating phase determines the downsamping phase
                            % of the equalizer output in function 'PolarizationDemux'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
for ktime = 1:nBlk
    X1 = [x1(:,ktime),x2(:,ktime)] .* h1;
    X2 = [x1(:,ktime),x2(:,ktime)] .* h2;
    Y1( ktime ) =  sum( sum( X1 ) ); % Calculating outputs
    Y2( ktime ) =  sum( sum( X2 ) );
    if stage == 1   % 1st stage CMA and update h1 only
        if (dontskip || mod(ktime-1,sps)==modk )
            incr1 = mu.*errorfuncma( Y1(ktime), R, err_id).*conj( xx(nindex,:) );
            h1 = h1 + incr1;
            h2(:,2) = ConjReverse( h1(:,1) );
            h2(:,1) = -1*ConjReverse( h1(:,2) );
        end
    else            % 2nd stage CMA and update h1h2
        if (dontskip || mod(ktime-1,sps)==modk)
            incr1 = mu.*errorfuncma( Y1(ktime), R, err_id).*conj( xx(nindex,:) );
            incr2 = mu.*errorfuncma( Y2(ktime), R, err_id).*conj( xx(nindex,:) );
            h1 = h1 + incr1;
            h2 = h2 + incr2;
        end
    end
    h = [h1(halftaps+1,:); h2(halftaps+1,:)];
    deth(ktime) = abs(det(h));
    switch method
        case 1  %% determination monitoring
            if abs(det(h)) < 0.01
                h1(:,1) = 0.5*( h1(:,1) + ConjReverse(h2(:,2)) );
                h2(:,2) = +ConjReverse( h1(:,1) );
                h1(:,2) = 0.5*( h1(:,2) - ConjReverse(h2(:,1)) );
                h2(:,1) = -ConjReverse( h1(:,2) );
            end
        case 2  %% alan method
            h1(:,1) = 0.5*( h1(:,1) + ConjReverse(h2(:,2)) );
            h2(:,2) = 0.5*( h2(:,2) + ConjReverse(h1(:,1)) );
            h1(:,2) = 0.5*( h1(:,2) - ConjReverse(h2(:,1)) );
            h2(:,1) = 0.5*( h2(:,1) - ConjReverse(h1(:,2)) );
    end
    if err_id == 1
        mse(ktime) = (abs(Y1(ktime)).^2 - R(1)).^2;
    else
        mse(ktime) = (real(Y1(ktime)).^2+imag(Y1(ktime)).^2-2*R(1)).^2;
    end
end
Y = [Y1 Y2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Err = errorfuncma( X, R, IDX)
switch IDX
    case 1   %%% Classical CMA
        Err = -X.*(abs(X).^2-R(1));
    case 2   %%% Modified CMA
        Err = -1.*complex(real(X).*(real(X).^2-R(1)),imag(X).*(imag(X).^2-R(2)));
    case 3   %%% Cascaded multi-modulus algorithm
        A1 = 0.5*(R(1)+R(2));
        A2 = 0.5*(R(3)-R(1));
        A3 = 0.5*(R(3)-R(2));
        e1 = abs(X)-A1; e2 = abs(e1)-A2; e3 = abs(e2)-A3;
        Err = -sign(X) * e3 * sign(e1) * sign(e2);
    case 4   %%% Radius directed
        Err = X.*min(R -abs(X).^2);
    case 5   %%% half-constellation MCMA
        window = windecision(real(X),2/3) || windecision(imag(X),2/3);
        Err = -1.*window.*complex(real(X).*(real(X).^2-R(1)),imag(X).*(imag(X).^2-R(2)));
    case 6   %%% Modified radius directed
        Err = complex(real(X).*min(R-real(X).^2),imag(X).*min(R-imag(X).^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ConjReverse(x)
y = conj(flipud(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = windecision(x, D)
y = 0.5*(1+sign(x.^2-D.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%