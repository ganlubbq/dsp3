function [Y h1 h2 mse deth] = CmaEqualizer2(xx, h1, h2, taps, mu, R, sps, err_id, stage, method)
%CMAADAPTIVEFILTER Polarization demultiplexing filter using CMA algorithm
%   [Y H1 H2] = CMAADAPTIVEFILTER(XX, H1, H2, TAPS, MU, R, SPS) applies a
%   matrix of adaptive filters whose initial state is written in H1 and H2
%   matrices. The filters coefficients are updated with the constant
%   modulus algorithm (CMA) proposed by Godard [1]. The parameters of this
%   algorithm are the number of taps TAPS, the radius R, the convergence
%   parameter MU and the number of samples per symbol SPS.
%   Normally there are four filters, let's call them H11, H12, H21, H22 and
%   the relation between the inputs and the outputs of this function is the
%   following:
%       Y1 = X1 ** H11 + X2 ** H12
%       Y2 = X1 ** H21 + X2 ** H22
%   where ** is the convolution, X1 and X2 are the two inputs signals
%   and Y1 and Y2 are the two output signals.
%   Usually the filters Hyx have several coefficients (taps).
%   In order to store all the filter a 3-dimensional matrix is needed. I
%   propose the following order for the dimensions: H( tap, y , x ) where H
%   has the following size: 5x2x2 if the filters have 5 taps each.
%   Since Matlab is very slow with operations with such a kind of matrix I
%   preferred to split H in two 2-dimensional matrices: h1 and h2 are the
%   matrix you could virtually obtain by doing:
%       H1 = squeeze( H( :, 1, : ) )
%       H2 = squeeze( H( :, 2, : ) )
%   So H11 = H1(:,1), H12 = H1(:,2), H21 = H2(:,1), H22 = H2(:,2)
%
%   [1] D. N. Godard, "Self-Recovering Equalization and Carrier Tracking in
%   Two-Dimensional Data Communication Systems," IEEE Trans. Commun., vol.
%   COM-28, no. 11, Nov. 1980
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
mse = zeros( 1 ,L );
deth = zeros( 1 ,L );
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