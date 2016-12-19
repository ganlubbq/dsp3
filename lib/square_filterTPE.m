% DESCRIPTION
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function [y,phEst] = square_filterTPE(x,symrate,sps,bs,decflag)
%square timing recovery, 4 samples/symbol in, 2 sps out
%

if nargin < 5
    decflag = 0;
end

IX = real(x(:,1));
QX = imag(x(:,1));
IY = real(x(:,2));
QY = imag(x(:,2));

SymbolRate = symrate;
SymbolPeriod = 1/SymbolRate;

if sps == 2
    
    Time2sps = (0: 1: length(IX)-1) * SymbolPeriod/2;
    Time4sps = (0:1:2*length(IX)-1) * SymbolPeriod/4;
    
    % InterpTechniqueDSF = 'interp1';
    InterpTechniqueDSF = 'interpft';
    
    % interpolation to get 4 samples per symbol using interp1
    if strcmp(InterpTechniqueDSF, 'interp1')
        % method = 'nearest';
        % method = 'linear';
        % method = 'spline';
        % method = 'pchip';
        method = 'cubic';
        IX_4sps = interp1(Time2sps, IX, Time4sps, method);
        QX_4sps = interp1(Time2sps, QX, Time4sps, method);
        IY_4sps = interp1(Time2sps, IY, Time4sps, method);
        QY_4sps = interp1(Time2sps, QY, Time4sps, method);
    end
    
    % interpolation to get 4 samples symbol bit using interpft
    if strcmp(InterpTechniqueDSF, 'interpft')
        IX_4sps = interpft(IX, 2*length(IX));
        QX_4sps = interpft(QX, 2*length(QX));
        IY_4sps = interpft(IY, 2*length(IY));
        QY_4sps = interpft(QY, 2*length(QY));
    end
end

if sps == 4
    IX_4sps = IX;
    QX_4sps = QX;
    IY_4sps = IY;
    QY_4sps = QY;
    Time4sps = (0:1:length(IX)-1) * SymbolPeriod/4;
end


NumBlk = fix(length(IX_4sps) / bs);

Xseq = abs(IX_4sps + 1i*QX_4sps).^2;
Yseq = abs(IY_4sps + 1i*QY_4sps).^2;

sumX = zeros( 1, NumBlk);
sumY = zeros( 1, NumBlk);

for k = 1 : NumBlk
    sum1 = 0;
    sum2 = 0;
    for m = 0 : bs-1
        sum1 = sum1 + Xseq((k-1)*bs + m+1) * exp(-1i*2*pi*m/4);
        sum2 = sum2 + Yseq((k-1)*bs + m+1) * exp(-1i*2*pi*m/4);
    end
    sumX(k) = sum1;
    sumY(k) = sum2;
end

% smooth the summation signal
NavgDSF = floor(length(sumX)/2);
filt_sumX = filter( ones(1,NavgDSF)/NavgDSF, 1, sumX);
filt_sumY = filter( ones(1,NavgDSF)/NavgDSF, 1, sumY);

% get the normalized delay
tauX_seq = -1/2/pi * unwrap(angle(filt_sumX));
tauY_seq = -1/2/pi * unwrap(angle(filt_sumY));

% determine the single value of delay
rX = sum(exp(1i*2*pi*tauX_seq));
rY = sum(exp(1i*2*pi*tauY_seq));
tauX_s = SymbolPeriod/2/pi * atan2(imag(rX), real(rX));
tauY_s = SymbolPeriod/2/pi * atan2(imag(rY), real(rY));
% tauX_m = SymbolPeriod * mean(tauX_seq);
% tauY_m = SymbolPeriod * mean(tauY_seq);

phEst = [ tauX_seq(:)*sps, tauY_seq(:)*sps];

% interpolate sampling based on estimated delays
TimeX_after = Time4sps + tauX_s;
TimeY_after = Time4sps + tauY_s;

method = 'linear';
% method = 'cubic';
IX_4spsRT = interp1(Time4sps, IX_4sps, TimeX_after, method);
QX_4spsRT = interp1(Time4sps, QX_4sps, TimeX_after, method);
IY_4spsRT = interp1(Time4sps, IY_4sps, TimeY_after, method);
QY_4spsRT = interp1(Time4sps, QY_4sps, TimeY_after, method);

I_X = downsample(IX_4spsRT, 2, 0);
Q_X = downsample(QX_4spsRT, 2, 0);
I_Y = downsample(IY_4spsRT, 2, 0);
Q_Y = downsample(QY_4spsRT, 2, 0);

% after timing recovery
signal_x_tim = I_X + 1i*Q_X;
signal_y_tim = I_Y + 1i*Q_Y;

y = [signal_x_tim(:),signal_y_tim(:)];

if decflag
    y = y(1:2:end,:);
end

return

