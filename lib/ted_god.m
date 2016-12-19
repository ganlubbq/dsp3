function y = ted_god(xx)

% g = 0.0001;
% L = length(xx);
% N = 256;            % FFT size
% bins = 64;          % Integral bins
% iter = floor(L/N);
% tau = zeros(iter,1);
% rr = reshape(real(xx(1:N*iter)),N,iter);
% % rr = reshape(xx(1:N*iter),N,iter);
% 
% fbaud = N/4;
% f_p = (fbaud-bins/2) : (fbaud+bins/2);
% f_baud = 3*N/4;
% f_n = (f_baud-bins/2) : (f_baud+bins/2);

% according to the Planchers principle, the production of rr is
% equivalent to that of RR.

% for k = 1:iter
%     RR  = fft(rr(:,k));
%     tau(k) = g*sum(imag(RR(f_p).*conj(RR(f_n))));
% end
% y = -mean(tau);



N = 2^nextpow2(length(xx));
yy = fft(xx,N);
% figure; semilogy(abs(yy).^2);


% sum over half of spectrum
xc = xcorr(yy);
% figure; semilogy(abs(xc).^2);
% 
% GCT
% y = -angle(xc(N/2)) / (2*pi);
% 
% % in Godard's original paper
y = -imag(xc(N/2)) / N;




% sum of single frequency bin at fbaud and -fbaud
% large tolerance to CD and PMD also large jitter
% xc = yy(N/4) * conj(yy(3*N/4));

% GCT
% y = -angle(xc) / (2*pi);

% % in Godard's original paper
% y = -imag(xc);



% % sum of few frequency bins at fbaud and -fbaud
% % large tolerance to CD and PMD also large jitter
% xc = yy(N/4-32:N/4+32) .* conj(yy(3*N/4-32:3*N/4+32));
% 
% % in Godard's original paper
% y = -imag(sum(xc));


