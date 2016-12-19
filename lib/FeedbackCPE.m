function [y, phi] = FeedbackCPE( signal, mn, mu, initial )
%FEEDFACKCPE Decision-directed feefback carrier phase recovery routine
% neglect the amplitude noise, data received rx = A*exp(j*phi1) and its
% decisiion dx = A*exp(j*phi0). The error signal err = rx*exp(-j*phierr)-dx
% = A*exp(j*phi0)*[exp(j(phi1-phi0-phierr))-1], where phierr is current
% estimated phase error and (phi1-phi0-phierr) would be the error of the
% estimated phase error. We take [rx*exp(-j*phierr)]*[conj(err)] =
% A^2*[1-exp(j*(phi1-phi0-phierr))], using imaginary part to approximate
% the error of phase error.
%
% CopyRight:Wang Dawei EIE PolyU   $Date:16/3/2010
N = length(signal);
y = zeros(size(signal));
phi= zeros(size(signal));
phi(1,:) = initial;
switch mn
    case 2
        for kk = 1:N
            zz = signal(kk,:) .* exp(-1j*phi(kk,:));
            y(kk,:) = zz;
            aa = sign(real(zz));
            err = zz - aa;
            phi(kk+1,:) = phi(kk,:) - mu*imag( zz.*conj(err) );
        end
    case 4
        for kk = 1:N
            zz = signal(kk,:) .* exp(-1j*phi(kk,:));
            y(kk,:) = zz;
            aa = sign(real(zz)) + 1j*sign(imag(zz));
            err = zz - aa;
            phi(kk+1,:) = phi(kk,:) - mu*imag( zz.*conj(err) );
        end
    case 16
        for kk = 1:N
            zz = signal(kk,:) .* exp(-1j*phi(kk,:));
            y(kk,:) = zz;
            aa = slicer16qam(zz);
            err = zz - aa;
            phi(kk+1,:) = phi(kk,:) - mu*imag( zz.*conj(err) );
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = slicer16qam(x)
x = x(:);
bound = 2;
a = [sign(real(x)), sign(real(x)-sign(real(x))*bound)];
a = a(:,2)+2*a(:,1);
b = [sign(imag(x)), sign(imag(x)-sign(imag(x))*bound)];
b = b(:,2)+2*b(:,1);
y = a + 1j*b;
y = y.';