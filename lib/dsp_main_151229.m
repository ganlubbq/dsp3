function [dspOut1, dspOut2, dspState] = dsp_main_151229(xi, xq, yi, yq, dspParam)
% Main testing funtion for offline processing. Using circular boundary
% condition for adaptive equalizer convergence. 
% 
% Project: digital communication recevier (DCR2015)
%
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: All units are international standard units
% 
% See Also: 
% 
% Copyright 2015 

dspState = [];
fprintf('+ processing %d input samples \n', 4 * length(xi));

% remove the dc componant
XI = xi - sum(xi) / length(xi);
XQ = xq - sum(xq) / length(xq);

YI = yi - sum(yi) / length(yi);
YQ = yq - sum(yq) / length(yq);

% unit the mean power
XI = XI / sqrt(sum(abs(XI).^2) / length(XI));
XQ = XQ / sqrt(sum(abs(XQ).^2) / length(XQ));

YI = YI / sqrt(sum(abs(YI).^2) / length(YI));
YQ = YQ / sqrt(sum(abs(YQ).^2) / length(YQ));

% remove the non-numerical
XI(isnan(XI)) = 0;
XQ(isnan(XQ)) = 0;

YI(isnan(YI)) = 0;
YQ(isnan(YQ)) = 0;

XI = XI(:);
YI = YI(:);
XQ = XQ(:);
YQ = YQ(:);

% the power spectral density
dataX = XI + 1i * XQ;
dataY = YI + 1i * YQ;

psdX = abs(fftshift(fft(dataX))) .^ 2;
psdY = abs(fftshift(fft(dataY))) .^ 2;

% ACR function is the ifft of PSD of signal
acrX = abs(ifft(psdX)) / length(psdX);
acrY = abs(ifft(psdY)) / length(psdY);

% figure; plot(acrX); hold; plot(acrY,'r'); hold off
% figure; plot(10*log10(psdX)); hold; plot(10*log10(psdY),'r'); hold off

if dspParam.showEye
    plotEyeDiagram(xi+1i*xq,dspParam.adcFs/dspParam.Rs*3,'e');
end

if dspParam.doDigitalLPF
    [xi, xq, yi, yq] = LowpassFilter(samplingRate,Pa_Bw,xi,xq,yi,yq);
end

REAL_DATA = [XI, XQ, YI, YQ];

if dspParam.showEye
    plotEyeDiagram(xi+1i*xq,dspParam.adcFs/dspParam.Rs*3,'e');
end

sps = 2;

% need to implement a realistic resampling algorithm here
fprintf('+ doing digital resampling \n');
[up, down] = rat(dspParam.Rs * sps / dspParam.adcFs);
if up == down
    % do nothing
elseif up == 1
	XI = XI(1 : down : end, 1);
	YI = YI(1 : down : end, 1);
	XQ = XQ(1 : down : end, 1);
	YQ = YQ(1 : down : end, 1);
else
    XI = resample(XI, up, down);
	YI = resample(YI, up, down);
	XQ = resample(XQ, up, down);
	YQ = resample(YQ, up, down);
end

fs = dspParam.Rs * sps;

rawX = XI + 1i * XQ;
rawY = YI + 1i * YQ;

if dspParam.doFrontEndComp
    fprintf('+ doing receiver front-end compensation \n');
    fecX = orthogonalization(rawX);
	fecY = orthogonalization(rawY);
else
    fecX = rawX;
	fecY = rawY;
end

nSample = length(fecX);

if dspParam.doCDE
    fprintf('+ doing chromatic dispersion estimation \n');
end

if dspParam.doCDC
    fprintf('+ doing chromatic dispersion compensation \n');
	[dspState.HCD] = calcDispResponse(nSample, fs, dspParam.lambda, dspParam.lambda0, dspParam.DL, dspParam.SL);
	cdcX = ifft(fft(fecX) .* conj(dspState.HCD));
	cdcY = ifft(fft(fecY) .* conj(dspState.HCD));
else
    cdcX = fecX;
	cdcY = fecY;
end

if dspParam.doTPE
    fprintf('+ doing timing recovery \n');
else
    tpeX = cdcX;
	tpeY = cdcY;
    dspState.TPE_PHASE = [];
end

if dspParam.doDownSampling
    fprintf('+ doing digital downsampling \n');
	tpeX = tpeX(1 : 2 : end);
	tpeY = tpeY(1 : 2 : end);
	sps = 1;
end

if dspParam.doMIMO
    fprintf('+ doing MIMO \n');
else
    cmaX = tpeX;
	cmaY = tpeY;
    dspState.CMA_MSE = [];
end

if dspParam.doFOE
    fprintf('+ doing frequency offset estimation \n');
else
    dspState.freqOffset = 0;
end

if dspParam.doCPE
    fprintf('+ doing carrier recovery | %s \n', dspParam.cpeAlgSelect);
	switch dspParam.cpeAlgSelect
		case 'BPS'
            dspState.CPE_PNx = estimateCarrierPhaseBPS(cmaX, dspParam.cpeBlockSize, dspParam.cpeBPSnTestPhase, dspParam.mn);
            dspState.CPE_PNy = estimateCarrierPhaseBPS(cmaY, dspParam.cpeBlockSize, dspParam.cpeBPSnTestPhase, dspParam.mn);
		case 'VVPE'
            dspState.CPE_PNx = estimateCarrierPhaseVV(cmaX, dspParam.cpeBlockSize, dspParam.vvpeAvgMode);
            dspState.CPE_PNy = estimateCarrierPhaseVV(cmaY, dspParam.cpeBlockSize, dspParam.vvpeAvgMode);
		case 3
			% TBA
		otherwise
			warning('unsupported cpe algorithm'); keyboard;
    end
    cpeX = cmaX .* exp(-1i * dspState.CPE_PNx);
    cpeY = cmaY .* exp(-1i * dspState.CPE_PNx);
else
    cpeX = cmaX;
	cpeY = cmaY;
    dspState.CPE_PNx = [];
	dspState.CPE_PNy = [];
end

if dspParam.doMLCPE
% 	cpeX = cpe_mle(cpeX,dspParam.cpeMlBlkSize,dspParam.cpeMlIter,dspParam.mn);
% 	cpeY = cpe_mle(cpeY,dspParam.cpeMlBlkSize,dspParam.cpeMlIter,dspParam.mn);
end

if dspParam.doLmsAfterCPE
    fprintf('+ doing LMS after carrier recovery \n');
    lmsX = lms_sng_eq(cpeX, dspParam.mn, sps, dspParam.lmsGainAfterCPE, dspParam.lmsTapsAfterCPE, dspParam.lmsIterAfterCPE);
    lmsY = lms_sng_eq(cpeY, dspParam.mn, sps, dspParam.lmsGainAfterCPE, dspParam.lmsTapsAfterCPE, dspParam.lmsIterAfterCPE);
else
	lmsX = cpeX;
	lmsY = cpeY;
end

dspOut1 = lmsX;
dspOut2 = lmsY;

fprintf('+ output %d samples \n', 2 * length(lmsX));

return
