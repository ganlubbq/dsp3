

function [y,H] = FrequencyDCF(x,Dz,fc,fs,ntaps,method,overlap)
%FREQUENCYDCF Frequency domain chromatic Dz compensation
%   Dz: the unit is ps/nm and is product of D and L
%   fc: center frequency
%   fs: sampling frequency
%   ntaps: number of FIR taps
%   method: 'overlap' or 'ideal'
%   overlap: ratio between overlap and one block
%
%   Example
%   
%   See also Fiber

%   Copyright2010 WANGDAWEI 16/3/2010

if nargin<7
    overlap = 0.5;
end

N = length(x);
Beta2L = -Dz * 299792458 / (2*pi*fc^2);

if Beta2L == 0
    y = x; H = []; return;
end

if strcmpi(method,'overlap')
    if ntaps > N
        ntaps = N;
    end
    % overlap length
    overlap = floor( ntaps * overlap);
    
    if mod(overlap,2)
        overlap = overlap+1;
    end
    
    % half the overlap, adjustable
    precursor = 0.5 * overlap;
    
    % transfer function of GVD
    H = DispCompFilter( Beta2L, fs, ntaps);
    H = H * ones(1,size(x,2));
    y = [];
    
    for ii = 1 : (ntaps-overlap) : (N-ntaps+1)
        % one block
        idx = ii : (ii+ntaps-1);
        % equalize
        z = ifft( fft(x(idx,:)) .* H );
        % cut off the precursor and postcursor
        if ii == 1
            z = z(1:(end-precursor),:);
        elseif ii == (N-ntaps+1)
            z = z((precursor+1):end,:);
        else
            z = z((precursor+1):(end-precursor),:);
        end
        y = [y;z];
    end
%     y = y(1:end-1,:);
elseif strcmpi(method,'fde')
    if ntaps > N
        ntaps = N;
    end
    % transfer function of GVD
    H = DispCompFilter( Beta2L, fs, ntaps);
    H = H * ones(1,size(x,2));
    y = [];
    for ii = 1 : ntaps : (N-ntaps+1)
        % one block
        idx = ii : (ii+ntaps-1);
        % equalize
        z = ifft( fft(x(idx,:)) .* H );
        y = [y;z];
    end
elseif strcmpi(method,'ideal')
    % transfer function of GVD
    H = DispCompFilter( Beta2L, fs, N);
    H = H * ones(1,size(x,2));
    % compensate
    y = ifft( fft(x) .* H );
end


function Hfilt = DispCompFilter(Beta2L, fs , N)
freq = (-fs/2 : fs/N : fs/2 * (N-2)/N);
argum = (2*pi*freq).^2 * Beta2L * 0.5;
H = exp ( 1j * argum);
Hfilt = ifftshift( H.');

% freq = linspace(0,fs,N);
% figure; 
% plotyy(freq,abs(fft(x)).^2, freq,angle(H(:,1)))
% xlabel('Frequency (Hz)')

% figure; plot(freq,a)
% xlabel( 'Frequency (Hz)' )
% ylabel( 'Phase' )
% title( 'GVD (all-pass) filter' )

