% Main testing funtion for offline processing. Using circular boundary
% condition for adaptive equalizer convergence.
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
% Copyright 2015 Default

function [dspOut1,dspOut2,dspState] = dsp_main_150927( xi, xq, yi, yq,dspParam)

dspState = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the dc componant
XI = xi - sum(xi)/length(xi);
XQ = xq - sum(xq)/length(xq);

YI = yi - sum(yi)/length(yi);
YQ = yq - sum(yq)/length(yq);

% unit the mean power
XI = XI/sqrt(sum(abs(XI).^2)/length(XI));
XQ = XQ/sqrt(sum(abs(XQ).^2)/length(XQ));

YI = YI/sqrt(sum(abs(YI).^2)/length(YI));
YQ = YQ/sqrt(sum(abs(YQ).^2)/length(YQ));

% remove the non-numerical
XI(isnan(XI)) = 0;
XQ(isnan(XQ)) = 0;

YI(isnan(YI)) = 0;
YQ(isnan(YQ)) = 0;

XI = XI(:);
YI = YI(:);
XQ = XQ(:);
YQ = YQ(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.showEye
    ShowEyediagram(xi,xq,symRate,samplingRate)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doDigitalLPF
    [xi xq yi yq] = LowpassFilter(samplingRate,Pa_Bw,xi,xq,yi,yq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL_DATA = [XI,XQ,YI,YQ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.showEye
    ShowEyediagram(xi,xq,symRate,samplingRate)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sps = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to implement a realistic resampling algorithm here
[up down] = rat( dspParam.Rs*sps / dspParam.adcFs);
if up == down
    XI = XI;
	YI = YI;
	XQ = XQ;
	YQ = YQ;
elseif up == 1
	XI = XI(1:down:end,1);
	YI = YI(1:down:end,1);
	XQ = XQ(1:down:end,1);
	YQ = YQ(1:down:end,1);
else
    XI = resample(XI, up, down);
	YI = resample(YI, up, down);
	XQ = resample(XQ, up, down);
	YQ = resample(YQ, up, down);
end

fs = dspParam.Rs*sps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPLEX_DATA = [RSP_TMP(:,1)+1j*RSP_TMP(:,2),RSP_TMP(:,3)+1j*RSP_TMP(:,4)];
rawX = XI + 1j*XQ;
rawY = YI + 1j*YQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doFrontEndComp
    fecX = DspAlg.Orthogonal(rawX);
	fecY = DspAlg.Orthogonal(rawY);
else
    fecX = rawX;
	fecY = rawY;
end

nSample = length(fecX);

if dspParam.doCDE
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doCDC
	[dspState.HCD] = calcDispResponse(nSample,fs,dspParam.lambda,dspParam.lambda0,dspParam.DL,dspParam.SL);
	cdcX = ifft(fft(fecX).*conj(dspState.HCD));
	cdcY = ifft(fft(fecY).*conj(dspState.HCD));
else
    cdcX = fecX;
	cdcY = fecY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doTPE
    set(handles.VISA_INFO_OUTPUT,'string',sprintf([str{1:2}])); pause(0.01);
    norFlag = 0;
    if strcmpi(TPE_estMeth,'lee') || strcmpi(TPE_estMeth,'none')
        [TPE_OUT,TPE_PHASE] = DspAlg.FeedforwardTPE(DCF_OUT,constPoint,sps, ...
            TPE_blk,TPE_bias,TPE_estMeth,TPE_intMeth,TPE_decFlag,norFlag);
    elseif strcmpi(TPE_estMeth,'gardner') || strcmpi(TPE_estMeth,'none')
        [TPE_OUT,TPE_PHASE] = DspAlg.FeedbackTPE(DCF_OUT,constPoint,sps, ...
            TPE_gain,TPE_estMeth,TPE_intMeth,TPE_decFlag,norFlag);
    end
else
    tpeX = cdcX;
	tpeY = cdcY;
    dspState.TPE_PHASE = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doDownSampling
	tpeX = tpeX(1:2:end);
	tpeY = tpeY(1:2:end);
	sps = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doMIMO
    polmux = 1;
    [CMA_OUT,dspState.CMA_MSE] = DspAlg.PolarizationDemux(TPE_OUT,constPoint,sps,...
        polmux,CMA_gain,CMA_taps,CMA_errID,CMA_iter,0);
	cmaX = CMA_OUT(:,1);
	cmaY = CMA_OUT(:,2);
else
    cmaX = tpeX;
	cmaY = tpeY;
    dspState.CMA_MSE = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doFOE
    [CMA_OUT,df] = DspAlg.FeedforwardFOC(CMA_OUT,symRate);
else
    df = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doCPE
	switch dspParam.cpeAlgSelect
		case 1
			[cpeX,dspState.CPE_PNx] = cpe_bps(cmaX,dspParam.cpeBlockSize,dspParam.cpeBPSnTestPhase,dspParam.mn);
			[cpeY,dspState.CPE_PNy] = cpe_bps(cmaY,dspParam.cpeBlockSize,dspParam.cpeBPSnTestPhase,dspParam.mn);
		case 2
			[cpeX,dspState.CPE_PNx] = cpe_vvpe(cmaX,dspParam.cpeBlockSize,dspParam.vvpeAvgMode);
			[cpeY,dspState.CPE_PNy] = cpe_vvpe(cmaY,dspParam.cpeBlockSize,dspParam.vvpeAvgMode);
		case 3
			% TBA
		otherwise
			error('unsupported cpe algorithm');
	end
else
    cpeX = cmaX;
	cpeY = cmaY;
    dspState.CPE_PNx = [];
	dspState.CPE_PNy = [];
end

if dspParam.doMLCPE
	cpeX = cpe_mle(cpeX,dspParam.cpeMlBlkSize,dspParam.cpeMlIter,dspParam.mn);
	cpeY = cpe_mle(cpeY,dspParam.cpeMlBlkSize,dspParam.cpeMlIter,dspParam.mn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dspParam.doLmsAfterCPE
    lmsX = lms_sng_eq( cpeX,dspParam.mn,sps,dspParam.lmsGainAfterCPE,dspParam.lmsTapsAfterCPE,dspParam.lmsIterAfterCPE);
    lmsY = lms_sng_eq( cpeY,dspParam.mn,sps,dspParam.lmsGainAfterCPE,dspParam.lmsTapsAfterCPE,dspParam.lmsIterAfterCPE);
else
	lmsX = cpeX;
	lmsY = cpeY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dspOut1 = lmsX;
dspOut2 = lmsY;

return


