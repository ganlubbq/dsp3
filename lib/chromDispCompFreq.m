function [y, H] = chromDispCompFreq(x, fs, lambda, lambda0, DL, SL, ntaps, method, overlap)
% Frequency domain chromatic dispersion compensation
%
% %   DL: accumulated dispersion ps/nm
% %   fc: center frequency
% %   fs: sampling frequency
% %   ntaps: FFT size
% %   method: 'overlap' or 'one-shot'
% %   overlap: overlap ratio in one block
%
%   Example
%   
%   See also 
if ~iscolumn(x)
    error('the first input has to be a column vector');
end

if nargin < 9
    overlap = 0.5;
end

N = length(x);
Beta2L = -DL * 299792458 / (2*pi*fc^2);

if Beta2L == 0
    y = x; H = []; return;
end

if ntaps > N
    ntaps = N;
end

if strcmpi(method,'overlap')
    % overlap length
    overlap = floor(ntaps * overlap);
    
    % make it even
    if mod(overlap,2)
        overlap = overlap + 1;
    end
    
    % half the overlap, adjustable
    precursor = 0.5 * overlap;
    
    % transfer function of GVD
    H = calcDispResponse(ntaps, fs, lambda, lambda0, DL, SL);
    y = [];
    
    for ii = 1 : (ntaps-overlap) : (N-ntaps+1)
        % one block
        idx = ii : (ii+ntaps-1);
        % equalize
        z = ifft(fft(x(idx,:)) .* H);
        % cut off the precursor and postcursor
        if ii == 1
            z = z(1:(end-precursor),:);
        elseif ii == (N-ntaps+1)
            z = z((precursor+1):end,:);
        else
            z = z((precursor+1):(end-precursor),:);
        end
        y = [y; z];
    end

elseif strcmpi(method,'non-overlap')
    % transfer function of GVD
    H = calcDispResponse(ntaps, fs, lambda, lambda0, DL, SL);
    y = [];
    for ii = 1 : ntaps : (N-ntaps+1)
        % one block
        idx = ii : (ii+ntaps-1);
        % equalize
        z = ifft(fft(x(idx,:)) .* H);
        y = [y; z];
    end
    
elseif strcmpi(method,'ideal')
    % transfer function of GVD
    H = calcDispResponse(length(x), fs, lambda, lambda0, DL, SL);
    % compensate
    y = ifft(fft(x) .* H);
end

return
