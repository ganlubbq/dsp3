function y = ted(px, mode)
%TED Timing error detector with various algorithms
%   Spectral-line methods for symbol synchronization based on nonlinear
%   operation of the input symbols include AVN, SLN, FLN, LOGN;
%   Spectral-correlation methods for symbol synchronization based on
%   spectrum autocorrelation includes DTP, SLN, LEE, GARDNER, GODARD;
N = length(px);
k = 1 : N;
ex = exp(-1i .* (k - 1) .* pi ./ 2);

if strcmpi(mode, 'AVN')
    s = sum(abs(px) .* ex.');
elseif strcmpi(mode, 'SLN')
    s = sum(abs(px).^2 .* ex.');
elseif strcmpi(mode, 'LOGN')
    f = log(1 + 10.*abs(px).^2);
    s = sum(f .* ex.');
elseif strcmpi(mode, 'FLN')
    s = sum(abs(px).^4 .* ex.');
end
y = -angle(s) / 2 / pi;

if strcmpi(alg, 'dtp')
    x = abs(px(1 : 2 : end - 1)) .^ 2;
    y = abs(px(2 : 2 : end)) .^ 2;
    % cosine
    y = mean((x - y) ./ sqrt(2));
    % sine
    % y = mean((y - x) ./ sqrt(2);
    
elseif strcmpi(alg, 'Gardner')
    px1 = px(1 : 2 : end - 2);
    px2 = px(2 : 2 : end - 1);
    px3 = px(3 : 2 : end);
    y = mean((real(px1) - real(px3)) .* real(px2) ...
        + (imag(px1) - imag(px3)) .* imag(px2));
    
elseif strcmpi(alg, 'Godard')
    nfft = 2^nextpow2(length(px));
    yy = fft(px, nfft);    
    % sum over half of spectrum
    xc = xcorr(yy);
    % in Godard's original paper
    y = -imag(xc(nfft / 2)) / nfft;
    
elseif strcmpi(alg, 'Lee')
    g = 1.414;
    n = 1 : length(px);
    ex1 = (-1).^(n - 1);
    sum_1 = sum(abs(px).^2 .* ex1.');
    
    ex2 = 1j * (-1).^(n-1);
    xh = px(2 : end);
    xx = px(1 : end - 1);
    ex2 = ex2(1 : end - 1);
    sum_2 = sum(real(conj(xx) .* xh) .* ex2.');
    
    sum_3 = g * sum_1 + sum_2;
    y = angle(sum_3) / 2 / pi;
    
elseif strcmpi(alg, 'SLN')
    N = length(px);
    x = interpft(px, 2*N);
    k = 1 : 2*N;
    ex = exp(-1j.*(k - 1) .* pi ./ 2);
    s = sum(abs(x).^2 .* ex.');    
    y = -angle(s) / 2 / pi;
    
else
    keyboard;
end


return
