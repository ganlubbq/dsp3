function H = calcFilterFreqResp(nsample, fs, order, bandwidth, type)
% Calculate the frequency response of digital filter with input order and
% 3dB bandwidth
%
% Example:
%   H = calcFilterFreqResp(nsample, fs, order, bandwidth, type)
%
% Input:
%       nsample     - number of samples
%       fs          - sampling frequency [Hz]
%       order       - filter order/roll-off factor
%       bandwidth   - filter 3dB bandwidth one-sided/baudrate [Hz]
%       type        - filter type
%
% Reference:
% [1] Martin Pfennigbauer et al "Choice of MUX/DEMUX filter
% characteristics for NRZ RZ and CSRZ..." JLT 2006
%
% Note: The PSD = H.^2
%
% See also: calcDispResponse, calcOptGaussFlt

if order < 0
    error('filter order has to be positive');
end

switch lower(type)
    case 'bessel'
        vAlpha = [1.0, 1.361654129, 1.755672389, 2.113917675, 2.427410702, 2.703395061];
        
        % frequency grid in Hz
        freqGrid = getFFTGrid(nsample, fs);
        
        % Laplace transformation
        s = 1i * freqGrid / bandwidth * vAlpha(order);
        
        % Order coeff
        n = order;
        for k = 0 : n
            b(k+1) = factorial(2*n - k) / (2^(n - k) * factorial(k) * factorial(n - k));
            D(k+1, :) = b(k+1) * s.^k;
        end
        
        % numerator/denominator
        H = D(1,:) ./ sum(D);
        
        % Remove group delay
        ndxFreq = find(freqGrid >= 0 & freqGrid < bandwidth / vAlpha(order));
        pp = polyfit(freqGrid(ndxFreq) / 1e9, unwrap(angle(H(ndxFreq))), 1);
        H = H .* exp(-1i * polyval(pp, freqGrid/1e9));
        
        % normalize
        H = H / abs(H(1));
        
    case 'gaussian'
        % frequency grid in Hz
        freqGrid = getFFTGrid(nsample, fs);
        
        % frequency response, amplitude spetrum
        H = exp(-0.5 * log(2) * (freqGrid / bandwidth) .^ (2 * order));
        
    case 'rc'
        freq = getFFTGridPos(nsample, fs);
        if ~mod(nsample, 2)
            freq = [freq, freq(end) + fs/nsample];
        end
        
        H = zeros(size(freq)) + eps;
        f_low = (1 - order) * 0.5 * bandwidth;
        f_high = (1 + order) * 0.5 * bandwidth;
        
        fndxl = find(freq <= f_low);
        fndxm = find(freq > f_low & freq <= f_high);
        
        H(fndxl) = ones(size(fndxl));
        H(fndxm) = 0.5 * (1 + cos(pi / order * (freq(fndxm) / bandwidth - 0.5 * (1 - order))));
        
        H = [H(1 : end - 1 + mod(nsample, 2)), conj(fliplr(H(2 : end)))];
        
    case 'rrc'
        H = sqrt(calcFilterFreqResp(nsample, fs, order, bandwidth, 'rc'));
        
    case 'nyquist'
        H = calcFilterFreqResp(nsample, fs, order, bandwidth, 'rc');
        
        freq = getFFTGrid(nsample, fs);
        f_high = (1 + order) * 0.5 * bandwidth;
        
        HSINC = ones(size(H));
        HSINC(abs(freq) <= f_high) = sinc(freq(abs(freq) <= f_high) ./ bandwidth);
        
        H = H ./ HSINC;
        
    case 'rect'
        freq = getFFTGrid(nsample, fs);
        H = zeros(size(freq)) + eps;
        H(abs(freq) <= bandwidth) = 1;
        
    otherwise
        keyboard;
end

H = H(:);

return
