function [y, phEst] = FeedforwardTPE(signal, mn, sps, szBlock, bias, estMeth, intMeth, decFlag, norFlag)
%FEEDFORWARDTPE Feedforward timing phase estimation and recovery. The input
%parameters are defined as follows:
%
%   szBlock: block size
%   estMeth: estimator method
%   intMeth: interpolation method
%   decFlag: decimate flag
%   norFlag: normalize flag
%
%   Example
%   
%   See also: FeedbackTPE

%   Copyright2012 wangdawei 16/6/2012

if nargin < 9
    norFlag = 1;
end
if nargin < 8
    decFlag = 1;
end
if nargin < 7
    intMeth = 'linear';
end
if nargin < 6
    estMeth = 'lee';
end
if nargin < 5
    bias = 1.0;
end

% if sps~=2
%     error('DSPALG::FF_TPE 2 sps required')
% end

if norFlag
    x = DspAlg.Normalize( signal, mn)/(sqrt(mn)-1);
else
    x = signal;
end

% get the size of input
mm = size(x,1);
kk = size(x,2);

% chose the estimator
switch estMeth
    case 'none'
        y = x;
        phEst = [];
    case 'lee'
        % make sure that the length of x can be diveded by block-size
        temp = mod(mm, szBlock);
        if temp
            x = [x; zeros(szBlock-temp,kk)];
        end
        for ii = 1:kk
            [y(:,ii) phEst(:,ii)]= Lee_PLL(x(:,ii),sps,bias,szBlock,intMeth);
        end
    case {'sln','fln','avn'}
        % make sure that the length of x can be diveded by block-size
        temp = mod(mm,szBlock);
        if temp
            x = [ x; zeros(szBlock-temp,kk)];
        end
        for ii = 1:kk
            [y(:,ii) phEst(:,ii)]= SLN_PLL(x(:,ii),sps,szBlock,intMeth,estMeth);
        end
    otherwise
        error('DSPALG::FF_TPE unsupported estimator method')
end

% decimate the output or not
if decFlag
    y = y( 1:sps:end, :);
end


function [symRtm phEst]= Lee_PLL(x,sps,bias,szBlk,intMeth)
%LEE_PLL This is NOT a phase lock loop
%

% test the sign of delay first
temp = angle(TED_lee(x,bias))/2/pi;

% skip the first sample if needed
if temp<0
    x = circshift(x,-1);
end

% parallelization
x = reshape(x,szBlk,[]);

% number of blocks
nBlock = size(x,2);

% output blk by blk
datInt = zeros( szBlk/2, 2*nBlock);

sumcomp = zeros(1,nBlock);

for kk = 1:nBlock
    oneBlk = x(:,kk);
    sumcomp(kk) = TED_lee(oneBlk,bias);
end

% smooth the complex signal
N_tmp = 128;
summation = filter(ones(1,N_tmp)/N_tmp,1,sumcomp);

% normalized delay
epsilon = angle(summation)/2/pi;

% normalized delay in samples
phEst = epsilon * sps;

for kk = 1:nBlock
    
    oneBlk = x(:,kk);
    
    % first nBlock columns store the correct sample
    kk1 = kk;
    kk2 = kk+nBlock;
    
    switch intMeth
        case 'none'
            % select without any interpolation
            if abs(epsilon(kk)) >= 0.25
                datInt(:,kk1) = oneBlk(2:2:end);
                datInt(:,kk2) = oneBlk(1:2:end);
            else
                datInt(:,kk1) = oneBlk(1:2:end);
                datInt(:,kk2) = oneBlk(2:2:end);
            end
        case 'linear'
            % each column represents one symbol
            buffer = reshape(oneBlk,2,[]);
            % interpolation coefficient
            intEff = phEst(kk) - floor(phEst(kk));
            datInt(:,kk1) = buffer(1,:) + intEff*(buffer(2,:)-buffer(1,:));
            datInt(:,kk2) = buffer(2,:) + intEff*(buffer(1,:)-buffer(2,:));
        case 'parabolic'
            % parabolic interpolation requires 3 sps at least
        otherwise
            error('DSPALG::FF_TPE_Lee_PLL unsupported interpolator method')
    end
end
% format output
sampSet1 = datInt(1:end/2);
sampSet2 = datInt(end/2+1:end);
sampSet = [sampSet1;sampSet2];
symRtm = sampSet(1:end);



function [y,phEst] = SLN_PLL(x,sps,bs,intMeth,estMeth)
%square timing recovery, 4 samples/symbol in, 2 sps out
%

IX = real(x);
QX = imag(x);

if sps == 2
    
    Time2sps = 1: 1: length(IX);
    Time4sps = (1: 1: 2*length(IX)) /2;
    
    % InterpTechniqueDSF = 'interp1';
    InterpTechniqueDSF = 'interpft';
    
    % interpolation to get 4 samples per symbol using interp1
    if strcmp(InterpTechniqueDSF, 'interp1')
        % method = 'linear';
        method = 'cubic';
        IX_4sps = interp1(Time2sps, IX, Time4sps, method);
        QX_4sps = interp1(Time2sps, QX, Time4sps, method);
    end
    
    % interpolation to get 4 samples per symbol using interpft
    if strcmp(InterpTechniqueDSF, 'interpft')
        IX_4sps = interpft(IX, 2*length(IX));
        QX_4sps = interpft(QX, 2*length(QX));
    end
end

if sps == 4
    IX_4sps = IX;
    QX_4sps = QX;
    Time4sps = 1:1:length(IX);
end

NumBlk = length(IX_4sps) / bs;

x4 = IX_4sps + 1i*QX_4sps;

% Xseq = abs(x4).^2;
x4blc = reshape(x4,bs,[]);

sumX = zeros( 1, NumBlk);

if strcmpi(estMeth,'sln')
    for k = 1 : NumBlk
        sumX(k) = ted_SLN(x4blc(:,k));
    end
end
if strcmpi(estMeth,'avn')
    for k = 1 : NumBlk
        sumX(k) = ted_AVN(x4blc(:,k));
    end
end
if strcmpi(estMeth,'fln')
    for k = 1 : NumBlk
        sumX(k) = ted_FLN(x4blc(:,k));
    end
end

% smooth signal
NavgDSF = 256;
filt_sumX = filter( ones(1,NavgDSF)/NavgDSF, 1, sumX);

epsilon = -1/2/pi * unwrap(angle(filt_sumX));
if strcmpi(estMeth,'fln')
    epsilon = epsilon - sign(epsilon) * 0.5;
end
phEst = epsilon(:) * sps;

phEst_tmp = repmat( phEst.',bs,1);
phEstExt = phEst_tmp(:);

% interpolate sampling based on estimated delays
TimeX_after = Time4sps + phEstExt.';

IX_4spsRT = interp1(Time4sps, IX_4sps, TimeX_after, intMeth);
QX_4spsRT = interp1(Time4sps, QX_4sps, TimeX_after, intMeth);

IX_4spsRT(isnan(IX_4spsRT)) = 0;
QX_4spsRT(isnan(QX_4spsRT)) = 0;

% I_X = downsample(IX_4spsRT, 2, 0);
% Q_X = downsample(QX_4spsRT, 2, 0);

% after timing recovery
signal_x_tim = IX_4spsRT + 1i*QX_4spsRT;

y = signal_x_tim(:);





function s = ted_FLN(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px).^4 .* ex.' );


function s = ted_SLN(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px).^2 .* ex.' );


function s = ted_AVN(px)
N = length(px);
k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);
s = sum( abs(px) .* ex.' );


function y = TED_lee(x,g)
L = length(x);
n = 1:L;
% cosine part
ex1 = (-1).^(n-1);          % ex1 = exp(-1j.*(ii-1).*pi);
sum_1 = sum( abs(x).^2 .* ex1.' );
% sine part
ex2 = 1j * (-1).^(n-1);     % ex2 = exp(-1j.*(ii-1.5).*pi);
xh = x(2:end);
xx = x(1:end-1);
ex2 = ex2(1:end-1);
sum_2 = sum( real(conj(xx).*xh) .* ex2.' );
% with biasing
y = g*sum_1 + sum_2;

