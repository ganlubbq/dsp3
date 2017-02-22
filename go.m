function vM = go(vSet)
%% Fiber-optics with DSP
%
% ----> datestr(datenum(now))
%
% Birth:   04-Sep-2015 14:10:02
% Rebuilt: 13-Mar-2016 06:43:24
% Rebuilt: 12 Nov-2016 15:52:00
% Rebuilt: 18-Feb-2017 17:19:42
%
% Data Streaming
%
% DSP mode: 0-Offline; 1-Real-time
%
% Refined: 18-Feb-2017 17:19:42

close all

%% Constant
LIGHT_SPEED             = 299792458;
BOLTZMAN                = 1.381e-23;
ELECTRON                = 1.602e-19;
FRAME_WINDOW            = 3;
TEMPERATURE             = 300;
PD_LOAD_RESISTANCE      = 5000; % TIA in PD
NOISE_REFERENCE_BAND    = 12.5e9; % for OSNR definition
DETECTION_MODE          = 'HOM'; % HOM or HET
EID                     = 'goErr';
CENTER_FREQUENCY        = 193.1e12;
CENTER_WAVELENGTH       = LIGHT_SPEED/CENTER_FREQUENCY;
LOG                     = 0;
FIG_TXPN                = 1;
FIG_TXLASER             = 2;
FIG_DRVEYE              = 3;
FIG_RECEIVED            = 4;
FIG_AFTERADC            = 5;

%% Switches
ctrlParam.doPilot       = 0;
ctrlParam.doNFC         = 0;
ctrlParam.doRZ          = 0;
ctrlParam.doNyquist     = 0;
ctrlParam.doRndPMD      = 0;
ctrlParam.doMzmComp     = 1;
ctrlParam.doCoherent    = 1;
ctrlParam.doDSP         = 1;
ctrlParam.doBalance     = 1;
ctrlParam.doPlot        = 0;
ctrlParam.verbose       = 1;

if nargin < 1 
    MAX_RUN_NUMBER          = 1000;
    HYBRID_90_PHASESHIFT    = 90; % degree
    ADC_SAMPLING_RATE       = 2; % samples per symbol
    DSP_MODE                = 0; % 0-offline; 1-real time
    DSO_MEMORY_LENGTH       = 64; % number of frames
    LASER_LINEWIDTH         = 500e3;
    OSNR                    = 150;
    baudrate                = 28e9;
    bitpersym               = 2;
    modFormat               = 'QPSK';
    freqOffset              = 1.0e9;
    psFiltType              = 'Nyquist'; % Nyquist Bessel Gaussian
else
    MAX_RUN_NUMBER          = vSet.nFrm;
    HYBRID_90_PHASESHIFT    = vSet.Hyd90; % degree
    ADC_SAMPLING_RATE       = vSet.ADCfs; % samples per symbol
    DSP_MODE                = vSet.DSPmode; % 0-offline; 1-real time
    DSO_MEMORY_LENGTH       = vSet.DSPmemLen; % number of frames
    LASER_LINEWIDTH         = vSet.linewidth;
    OSNR                    = vSet.osnr;
    baudrate                = vSet.buffer.rxPhaseNoise;
    bitpersym               = vSet.bitpersym;
    modFormat               = vSet.modFormat;
end

%[simulation]
samplingFs          = 8 * baudrate;
timewindow          = 512 / baudrate;
symbolsPerFrame     = timewindow * baudrate;
samplesPerSym       = samplingFs / baudrate;
samplesPerFrame     = symbolsPerFrame * samplesPerSym;
vctFreqPerFrm       = getFFTGrid(samplesPerFrame, samplingFs);
numSamples          = FRAME_WINDOW * samplesPerFrame;
timeVector          = (0:numSamples-1)' / samplingFs;
freqVector          = getFFTGrid(numSamples, samplingFs);
ALPHABET_SIZE       = 2^bitpersym;
vM.StartTime        = datestr(now);
timeVectorAbs       = timeVector; % initialize absolute time

if LOG
    logFile = fopen('log','a');
    fprintf(logFile,'-------------------------------------');
    fprintf(logFile,'\n');
    fprintf(logFile,'Time now is %s\n',datestr(now));
    if exist('vSet','var')
        fprintf(logFile,'Caller is %s\n', vSet.caller);
    end
else
    logFile = 1;
end


%% Project setting
Txsamplerate        = samplingFs;

txLaser.centerFreq = CENTER_FREQUENCY;
txLaser.centerLambda = LIGHT_SPEED / txLaser.centerFreq;
txLaser.linewidth = LASER_LINEWIDTH;
txLaser.phaseNoiseVar = 2*pi * txLaser.linewidth / samplingFs;
txLaser.azimuth = 0;
txLaser.ellipticity = 0;
txLaser.power = 1e-3; % [W]
txLaser.RIN = -130; % dBc/Hz

modulator.Vpi = 3; % [V]
modulator.extRatio = 35; % [dB]
modulator.efficiency            = 0.55; % modulation efficiency

transmtter.orderLPF = 4;
transmtter.bandwidth = 0.75 * baudrate; % [Hz]
transmtter.pilotX = [];
transmtter.pilotY = [];

fiber.n2 = 2.6e-20;
fiber.coreArea = 80e-12; % [m^2]
% fiber.initAngle = 0; % [degree]
fiber.lossFast = 0.2; % [dB/km]
fiber.lossSlow = 0.2;
fiber.spanLength = 80e3; % [m]
fiber.spanNum = 2;
chn_stepLength = 1e3; % [m]
chn_corrLength = 100; % [m]
chn_dispD = 17e-6; % [s/m]
chn_dispS = 0.08e3; % [s/m^2]
chn_fullPMD = 0;
chn_PMDparam = 0.5e-12/31.623;
chn_ASEseed = 0;
optFilterOrder = 4;
optFilterBw = 40e9;
smfPolarIniAngle = 0;
smfPolarRotSpeed = 30e3; % rad/s
lnk_DL = chn_dispD*fiber.spanLength*fiber.spanNum;
lnk_SL = chn_dispS*fiber.spanLength*fiber.spanNum;

Rxcenterfreq        = CENTER_FREQUENCY;
RxcenterWave        = LIGHT_SPEED/Rxcenterfreq;
RxlaserSeed         = 0;
Rxlinewidth         = LASER_LINEWIDTH;
Rxsamplerate        = samplingFs;
Rxpnvar             = 2*pi*Rxlinewidth/Rxsamplerate;
rxLpfOrder          = 5;
rxLpfBw             = 0.75*baudrate;
RxLaserAzi          = 45;
RxLaserEll          = 0;
rxLaserPow          = 10e-3;
rxLaserRIN          = -130; % dB/Hz
rxPbcPowSpltRatio   = 0.5; % power split ratio of PBC
rxPbcPhaseRetard    = 0; % phase retardation of PBC
rxPdR               = 1.0;
rxPdDark            = 10E-9;

pd_thmvar   = 4 * BOLTZMAN * TEMPERATURE / PD_LOAD_RESISTANCE * (0.5 * samplingFs);




% preparing filter responses
txPulseShapeFilter.RollOffFactor = 0.35;
txPulseShapeFilter.freqRespRC = calcRCFreqResponse(numSamples, samplingFs, baudrate, txPulseShapeFilter.RollOffFactor, 0);
txPulseShapeFilter.freqRespRRC = calcRCFreqResponse(numSamples, samplingFs, baudrate, txPulseShapeFilter.RollOffFactor, 1);
txPulseShapeFilter.freqRespBessel = calcBesselResponse(numSamples, samplingFs, transmtter.orderLPF, transmtter.bandwidth);
txPulseShapeFilter.freqRespNyquist = calcNyquistFiltResponse(numSamples, samplingFs, baudrate, 0.1, 0);
txPulseShapeFilter.freqRespGaussian = calcGaussFlt(numSamples, samplingFs, 0, transmtter.orderLPF, transmtter.bandwidth);

rxPulseShapeFilter.freqRespBessel = calcBesselResponse(numSamples, samplingFs, rxLpfOrder, rxLpfBw);
rxPulseShapeFilter.freqRespNyquist = calcNyquistFiltResponse(numSamples, samplingFs, baudrate, 0.1, 0);
rxPulseShapeFilter.freqRespGaussian = calcGaussFlt(numSamples, samplingFs, 0, rxLpfOrder, rxLpfBw);

optGaussianFilter = calcOptGaussFlt(numSamples, samplingFs, 0, optFilterOrder, optFilterBw);


% initialize buffer
buffer.txBitsX = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.txBitsY = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.rxBitsX = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.rxBitsY = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.txPhaseNoise = genLaserPhaseNoise(numSamples, txLaser.phaseNoiseVar, 0);
buffer.rxPhaseNoise = genLaserPhaseNoise(numSamples, Rxpnvar, 0);
buffer.txLaserRIN   = zeros(1,numSamples);
buffer.rxLaserRIN   = zeros(1,numSamples);
buffer.channelNoiseX  = zeros(1,numSamples);
buffer.channelNoiseY  = zeros(1,numSamples);
buffer.thermNsIx    = zeros(1,numSamples);
buffer.thermNsQx    = zeros(1,numSamples);
buffer.thermNsIy    = zeros(1,numSamples);
buffer.thermNsQy    = zeros(1,numSamples);
buffer.shotNsIx     = zeros(1,numSamples);
buffer.shotNsQx     = zeros(1,numSamples);
buffer.shotNsIy     = zeros(1,numSamples);
buffer.shotNsQy     = zeros(1,numSamples);


% clear memory
memory1				= [];
memory2				= [];
memory3				= [];
memory4				= [];

dspMemCount		= 0;

% Controller
sysParam.addCD              = 0;
sysParam.addPMD             = 0;
sysParam.addLaserRIN        = 0;
sysParam.addLaserPN         = 0;
sysParam.addPolarRot        = 0;
sysParam.addThermNoise      = 0;
sysParam.addShotNoise       = 0;
sysParam.addAdcClockJitter  = 0;
sysParam.addFreqOffset      = 0;
sysParam.doBalanced         = 1;
sysParam.addFreqOffset      = 0;

% DSP
dspParam.Rs					= 28e9;
dspParam.mn					= 4;
dspParam.adcFs				= 56e9;
dspParam.showEye			= 0;
dspParam.doFrontEndComp		= 0;
dspParam.doDigitalLPF		= 0;
dspParam.doCDC				= 0;
dspParam.doCDE				= 0;
dspParam.doFramer			= 0;
dspParam.doTraining			= 0;
dspParam.doMIMO				= 0;
dspParam.doTPE				= 0;
dspParam.doCPE				= 1;
dspParam.doFOE				= 0;
dspParam.doFEC				= 0;
dspParam.doDownSampling = 1;
dspParam.DL = lnk_DL;
dspParam.SL = lnk_SL;
dspParam.lambda = CENTER_WAVELENGTH;
dspParam.lambda0 = CENTER_WAVELENGTH;
dspParam.lpfBW = 0;
dspParam.mimoStep = 1e-3;
dspParam.mimoTapNum = 7; %7 or 13
dspParam.mimoErrId = 0;
dspParam.mimoIteration = 10;
dspParam.doLmsAfterCPE = 0;
dspParam.lmsGainAfterCPE = 1e-3;
dspParam.lmsTapsAfterCPE = 125;
dspParam.lmsIterAfterCPE = 4;
dspParam.cpeBlockSize = 16;
dspParam.pllMu1 = 1e-3;
dspParam.pllMu2 = 1e-6;
dspParam.doMLCPE = 2;
dspParam.cpeMlBlkSize = 32;
dspParam.cpeMlIter = 1;
dspParam.vvpeAvgMode = 1; %0 for block 1 for sliding window
dspParam.cpeBPSnTestPhase = 10;
dspParam.cpeAlgSelect = 1; %1 for bps 2 for vvpe
refBitsOfflineX		= [];
refBitsOfflineY		= [];


%% Start main loop
for RUN = 1:MAX_RUN_NUMBER

% Binary source for a new frame
bitsX = randi([0 1],bitpersym,symbolsPerFrame);
bitsY = randi([0 1],bitpersym,symbolsPerFrame);

% fifo
[buffer.txBitsX] = fifoBuffer(buffer.txBitsX,bitsX);
[buffer.txBitsY] = fifoBuffer(buffer.txBitsY,bitsY);

% offline mode; take the center frame as the reference
if DSP_MODE == 0
	refBitsOfflineX = [refBitsOfflineX buffer.txBitsX(:,(1:symbolsPerFrame)+symbolsPerFrame)];
	refBitsOfflineY = [refBitsOfflineY buffer.txBitsY(:,(1:symbolsPerFrame)+symbolsPerFrame)];
elseif DSP_MODE == 1
	refBitsX = buffer.txBitsX(:,1:symbolsPerFrame);
	refBitsY = buffer.txBitsY(:,1:symbolsPerFrame);
end


%% Bit mapping
% convert bits to syms
txBaudX = symbolizerGrayQam(buffer.txBitsX);
txBaudY = symbolizerGrayQam(buffer.txBitsY);

% polarizationAnalyzer(txBaudX.',txBaudY.','o-');


%% TX laser pol x/y 

% generate new phase noise for the new frame
tmpPN = genLaserPhaseNoise(samplesPerFrame, txLaser.phaseNoiseVar, buffer.txPhaseNoise(end));

% fifo
buffer.txPhaseNoise = fifoBuffer(buffer.txPhaseNoise, tmpPN);

% buffer laser RIN is really defined as one-sided power variance

% generate new rin for new frame
tmpRIN = genWGN(1,samplesPerFrame, idbw(txLaser.RIN) * (txLaser.power^2) * (0.5*samplingFs), 'linear', 'real');

% fifo
buffer.txLaserRIN = fifoBuffer(buffer.txLaserRIN, tmpRIN);

if ctrlParam.doPlot
    pltIndex = 1:500;
    figure(FIG_TXPN); plot(tmpPN(pltIndex)); grid on; title('Tx Phase Noise');
end

% laser on
if sysParam.addLaserRIN
    if sysParam.addLaserPN
        txLaser.wfm = sqrt(txLaser.power + buffer.txLaserRIN) .* exp(1j * buffer.txPhaseNoise);
    else
        txLaser.wfm = sqrt(txLaser.power + buffer.txLaserRIN);
    end
else
    if sysParam.addLaserPN
        txLaser.wfm = sqrt(txLaser.power) .*  exp(1j * buffer.txPhaseNoise);
    else
        txLaser.wfm = sqrt(txLaser.power) * ones(1,numSamples);
    end
end

% if ctrlParam.doPlot
%     figure(FIG_TXLASER); plot(abs(txLaser(pltIndex).^2)); grid on; title('Tx laser power waveform');
% end


%% Driver
% normalize
txBaudRealX = real(txBaudX) / (sqrt(ALPHABET_SIZE)-1);
txBaudImagX = imag(txBaudX) / (sqrt(ALPHABET_SIZE)-1);
txBaudRealY = real(txBaudY) / (sqrt(ALPHABET_SIZE)-1);
txBaudImagY = imag(txBaudY) / (sqrt(ALPHABET_SIZE)-1);

% pre-distortion
if ctrlParam.doMzmComp
    txDrvIx = asin(txBaudRealX) * modulator.Vpi /pi;
    txDrvQx = asin(txBaudImagX) * modulator.Vpi /pi;
    txDrvIy = asin(txBaudRealY) * modulator.Vpi /pi;
    txDrvQy = asin(txBaudImagY) * modulator.Vpi /pi;
else
    txDrvIx = (txBaudRealX * modulator.efficiency) * modulator.Vpi /pi;
    txDrvQx = (txBaudImagX * modulator.efficiency) * modulator.Vpi /pi;
    txDrvIy = (txBaudRealY * modulator.efficiency) * modulator.Vpi /pi;
    txDrvQy = (txBaudImagY * modulator.efficiency) * modulator.Vpi /pi;
end


%% MZM non-linear pre-comp
% *add MZM nonlinearity pre-comp. here...*


%% Pulse shaping
% DAC - simple oversampling by inserting zeros
txDrvIxUps = upSampInsertZeros(txDrvIx, samplesPerSym);
txDrvQxUps = upSampInsertZeros(txDrvQx, samplesPerSym);

txDrvIyUps = upSampInsertZeros(txDrvIy, samplesPerSym);
txDrvQyUps = upSampInsertZeros(txDrvQy, samplesPerSym);


% bessel filtering
txDrvIxWfm = real(ifft(fft(txDrvIxUps) .* txPulseShapeFilter.freqRespRC));
txDrvQxwfm = real(ifft(fft(txDrvQxUps) .* txPulseShapeFilter.freqRespRC));
txDrvIyWfm = real(ifft(fft(txDrvIyUps) .* txPulseShapeFilter.freqRespRC));
txDrvQyWfm = real(ifft(fft(txDrvQyUps) .* txPulseShapeFilter.freqRespRC));

if ctrlParam.doPlot
    if exist('FIG_DRVEYE','var')
        plotEyeDiagram(FIG_DRVEYE,txDrvIxWfm, baudrate, samplingFs, 'electrical');
    else
        FIG_DRVEYE = plotEyeDiagram([],txDrvIxWfm, baudrate, samplingFs, 'electrical');
    end
    title('Electrical dirver signal');
end


% keyboard;
% % cr
% dat_drv_xi = clockRecIdeal(dat_drv_xi,txDrvIxWfm);
% dat_drv_xq = clockRecIdeal(dat_drv_xq,txDrvQxWfm);
% dat_drv_yi = clockRecIdeal(dat_drv_yi,txDrvIyWfm);
% dat_drv_yq = clockRecIdeal(dat_drv_yq,txDrvQyWfm);


%% MZM

% the bias point
V1 = - modulator.Vpi/2;
V2 = - modulator.Vpi/2;
V3 = + modulator.Vpi/2;

txOptSigX = oeModIqNested(txLaser.wfm, txDrvIxWfm, txDrvQxwfm, modulator.extRatio, modulator.Vpi, V1, V2, V3);
txOptSigY = oeModIqNested(txLaser.wfm, txDrvIyWfm, txDrvQyWfm, modulator.extRatio, modulator.Vpi, V1, V2, V3);


if ctrlParam.doPlot
    if exist('FIG_TXEYE','var')
        plotEyeDiagram(FIG_TXEYE,txOptSigX, baudrate, samplingFs, 'optical');
    else
        FIG_TXEYE = plotEyeDiagram([],txOptSigX, baudrate, samplingFs, 'optical');
    end
    title('Optical eye after modulator');
end


%% OSNR
% here goes an OSNR emulator



%% Fiber channel
% simulating the real fiber
if ~ctrlParam.doRndPMD % if random birefrigence is switched off, use simple model
    
    % convert osnr to snr per sample
    SNR = osnr2snr(OSNR, baudrate, samplesPerSym, 'complex');
    
    % calc noise power, should we use the power of optimal samples only ??
    sigPowXdb   = dbw(calcrms(txOptSigX)^2);
    sigPowYdb   = dbw(calcrms(txOptSigY)^2);
    noisePowerx = idbw(sigPowXdb - SNR);
    noisePowery = idbw(sigPowYdb - SNR);
    
    % generate new noise for new frame 
    channelNoiseX = genWGN(1,samplesPerFrame,noisePowerx,'linear','complex');
    channelNoiseY = genWGN(1,samplesPerFrame,noisePowery,'linear','complex');
    
    % fifo
    buffer.channelNoiseX = fifoBuffer(buffer.channelNoiseX, channelNoiseX);
    buffer.channelNoiseY = fifoBuffer(buffer.channelNoiseY, channelNoiseY);
    
    % add noise
    txOptSigX = txOptSigX + buffer.channelNoiseX(:);
    txOptSigY = txOptSigY + buffer.channelNoiseY(:);
    
    % fiber
    lk_rotjump = mod(smfPolarRotSpeed * timeVectorAbs, 2*pi);
    lk_theta = smfPolarIniAngle + lk_rotjump;
%     smfPolarIniAngle = lk_theta(numSamples);
    
    if sysParam.addPolarRot % rotate device only if polarization rotation is on
        tmpOptSigX = txOptSigX.*cos(lk_theta) - txOptSigY.*sin(lk_theta);
        tmpOptSigY = txOptSigX.*sin(lk_theta) + txOptSigY.*cos(lk_theta);
    else
        tmpOptSigX = txOptSigX;
        tmpOptSigY = txOptSigY;
    end
    
    if sysParam.addCD        
        [HCD] = calcDispResponse(numSamples, samplingFs, centerWave, centerWave, lnk_DL, lnk_SL);
        tmpOptSigX = ifft(fft(tmpOptSigX) .* HCD);
        tmpOptSigY = ifft(fft(tmpOptSigY) .* HCD);
    end
    
    if sysParam.addPMD
        tmpOptSigX = ifft(fft(tmpOptSigX) .* exp(-1j*2*pi*freqVector*(+lk_DGD/2)).*exp(-1j*2*pi*uv_fc*(+lk_DGD/2)));
        tmpOptSigY = ifft(fft(tmpOptSigY) .* exp(-1j*2*pi*freqVector*(-lk_DGD/2)).*exp(-1j*2*pi*uv_fc*(-lk_DGD/2)));
    end
    
    if sysParam.addPolarRot % rotate back only if polarization rotation is on 
        tmpOptSigX = tmpOptSigX.*cos(-lk_theta) - tmpOptSigY.*sin(-lk_theta);
        tmpOptSigY = tmpOptSigX.*sin(-lk_theta) + tmpOptSigY.*cos(-lk_theta);
    end
    
    if sysParam.addPolarRot % add endless rotation
        tmpOptSigX = tmpOptSigX.*cos(lk_theta) - tmpOptSigY.*sin(lk_theta);
        tmpOptSigY = tmpOptSigX.*sin(lk_theta) + tmpOptSigY.*cos(lk_theta);
        % update absolute time
        timeVectorAbs = ((0:numSamples-1)' + samplesPerFrame) / samplingFs;
        lk_rotjump = mod(smfPolarRotSpeed * timeVectorAbs, 2*pi);
    end
    
else % if random birefrigence is on, use SSF model
    
    % set initial OSNR
    tmpOptSigX = txOptSigX;
    tmpOptSigY = txOptSigY;
    
    % fiber input power control
end

rxOptSigX = tmpOptSigX;
rxOptSigY = tmpOptSigY;

if ctrlParam.doPlot
    % plot spectrum before and after fiber
end


%% RX laser
% polarization control
loJonesVector = calcJonesVector(RxLaserAzi, RxLaserEll);

% generate new phase noise for new frame
tmpPN = genLaserPhaseNoise(samplesPerFrame, Rxpnvar, buffer.rxPhaseNoise(end));

% fifo
buffer.rxPhaseNoise = fifoBuffer(buffer.rxPhaseNoise, tmpPN);

% buffer laser rin
tmpRIN = genWGN(1,samplesPerFrame,idbw(rxLaserRIN)*(rxLaserPow^2)*(samplingFs/2),'linear','real');
buffer.rxLaserRIN = fifoBuffer(buffer.rxLaserRIN, tmpRIN);

% laser on
if sysParam.addLaserRIN
    if sysParam.addLaserPN
        rxLaser = sqrt(rxLaserPow + buffer.rxLaserRIN) .* exp(1j*buffer.rxPhaseNoise);
    else
        rxLaser = sqrt(rxLaserPow + buffer.rxLaserRIN);
    end
else
    if sysParam.addLaserPN
        rxLaser = sqrt(rxLaserPow) .* exp(1j*buffer.rxPhaseNoise);
    else
        rxLaser = sqrt(rxLaserPow) * ones(1,length(buffer.rxLaserRIN));
    end
end

% control the polarization
rxLaser = loJonesVector * sqrt(abs(rxLaser).^2) .* [sign(rxLaser);sign(rxLaser)];

% PBS / PBC
loLaserPx = sqrt(rxPbcPowSpltRatio) * exp(-1j*rxPbcPhaseRetard) * rxLaser(1,:).';
loLaserPy = sqrt(1-rxPbcPowSpltRatio) * rxLaser(2,:).';


%% Hybrid
hybrid90 = deg2rad(HYBRID_90_PHASESHIFT);

switch DETECTION_MODE
    case 'HOM'
        % freq_offset     = Rxcenterfreq - txLaser.centerFreq;
        if sysParam.addFreqOffset
            xRealP = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .*  loLaserPx;
            xRealN = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .* -loLaserPx;
            xImagP = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .*  exp(1j*hybrid90) .*loLaserPx;
            xImagN = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .* -exp(1j*hybrid90) .*loLaserPx;
            yRealP = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .*  loLaserPy;
            yRealN = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .* -loLaserPy;
            yImagP = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .*  exp(1j*hybrid90) .*loLaserPy;
            yImagN = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .* -exp(1j*hybrid90) .*loLaserPy;
        else
            xRealP = rxOptSigX + loLaserPx;
            xRealN = rxOptSigX - loLaserPx;
            xImagP = rxOptSigX + exp(1j*hybrid90) .*loLaserPx;
            xImagN = rxOptSigX - exp(1j*hybrid90) .*loLaserPx;
            yRealP = rxOptSigY + loLaserPy;
            yRealN = rxOptSigY - loLaserPy;
            yImagP = rxOptSigY + exp(1j*hybrid90) .*loLaserPy;
            yImagN = rxOptSigY - exp(1j*hybrid90) .*loLaserPy;
        end
    case 'HET'
        if sysParam.addFreqOffset
            xHetP = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .*  loLaserPx;
            xHetN = rxOptSigX + exp(1j*2*pi*freqOffset*timeVector) .* -loLaserPx;
            yHetP = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .*  loLaserPy;
            yHetN = rxOptSigY + exp(1j*2*pi*freqOffset*timeVector) .* -loLaserPy;
        else
            xHetP = rxOptSigX + loLaserPx;
            xHetN = rxOptSigX - loLaserPx;
            yHetP = rxOptSigY + loLaserPy;
            yHetN = rxOptSigY - loLaserPy;
        end
    otherwise
        error('detection mode not supported !!');
end


%% Photo detector
% dark current
IpdDark = rxPdDark;

switch DETECTION_MODE
    case 'HOM'
        % square law detection
        IpdXrealP = rxPdR .* abs(xRealP).^2;
        IpdXrealN = rxPdR .* abs(xRealN).^2;
        IpdXimagP = rxPdR .* abs(xImagP).^2;
        IpdXimagN = rxPdR .* abs(xImagN).^2;
        IpdYrealP = rxPdR .* abs(yRealP).^2;
        IpdYrealN = rxPdR .* abs(yRealN).^2;
        IpdYimagP = rxPdR .* abs(yImagP).^2;
        IpdYimagN = rxPdR .* abs(yImagN).^2;
        
        if sysParam.addThermNoise
            IpdThermal1 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal2 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal3 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal4 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal5 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal6 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal7 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
            IpdThermal8 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2), pd_thmvar,'linear','real');
        else
            IpdThermal1 = 0;
            IpdThermal2 = 0;
            IpdThermal3 = 0;
            IpdThermal4 = 0;
            IpdThermal5 = 0;
            IpdThermal6 = 0;
            IpdThermal7 = 0;
            IpdThermal8 = 0;
        end
        
        if sysParam.addShotNoise
            pd_shtvar1   = 2*ELECTRON*(calcrms(xRealP)^2)*(samplingFs/2); % assuming responsitivity is 1
            pd_shtvar2   = 2*ELECTRON*(calcrms(xRealN)^2)*(samplingFs/2);
            pd_shtvar3   = 2*ELECTRON*(calcrms(xImagP)^2)*(samplingFs/2);
            pd_shtvar4   = 2*ELECTRON*(calcrms(xImagN)^2)*(samplingFs/2);
            pd_shtvar5   = 2*ELECTRON*(calcrms(yRealP)^2)*(samplingFs/2);
            pd_shtvar6   = 2*ELECTRON*(calcrms(yRealN)^2)*(samplingFs/2);
            pd_shtvar7   = 2*ELECTRON*(calcrms(yImagP)^2)*(samplingFs/2);
            pd_shtvar8   = 2*ELECTRON*(calcrms(yImagN)^2)*(samplingFs/2);
            IpdShot1 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar1,'linear','real');
            IpdShot2 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar2,'linear','real');
            IpdShot3 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar3,'linear','real');
            IpdShot4 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar4,'linear','real');
            IpdShot5 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar5,'linear','real');
            IpdShot6 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar6,'linear','real');
            IpdShot7 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar7,'linear','real');
            IpdShot8 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar8,'linear','real');
        else
            IpdShot1 = 0;
            IpdShot2 = 0;
            IpdShot3 = 0;
            IpdShot4 = 0;
            IpdShot5 = 0;
            IpdShot6 = 0;
            IpdShot7 = 0;
            IpdShot8 = 0;
        end
        
        V1 = IpdXrealP + IpdDark + IpdThermal1 + IpdShot1;
        V2 = IpdXrealN + IpdDark + IpdThermal2 + IpdShot2;
        V3 = IpdXimagP + IpdDark + IpdThermal3 + IpdShot3;
        V4 = IpdXimagN + IpdDark + IpdThermal4 + IpdShot4;
        V5 = IpdYrealP + IpdDark + IpdThermal5 + IpdShot5;
        V6 = IpdYrealN + IpdDark + IpdThermal6 + IpdShot6;
        V7 = IpdYimagP + IpdDark + IpdThermal7 + IpdShot7;
        V8 = IpdYimagN + IpdDark + IpdThermal8 + IpdShot8;
        
        if sysParam.doBalanced % balanced detection
            pd_xi = V1 - V2;
            pd_xq = V3 - V4;
            pd_yi = V5 - V6;
            pd_yq = V7 - V8;
        else
            pd_xi = V1;
            pd_xq = V3;
            pd_yi = V5;
            pd_yq = V7;
        end
        % low-pass filtering
        pd_xi = real(ifft(fft(pd_xi).* rxPulseShapeFilter.freqRespBessel));
        pd_xq = real(ifft(fft(pd_xq).* rxPulseShapeFilter.freqRespBessel));
        pd_yi = real(ifft(fft(pd_yi).* rxPulseShapeFilter.freqRespBessel));
        pd_yq = real(ifft(fft(pd_yq).* rxPulseShapeFilter.freqRespBessel));
        
    case 'HET'
        % square law detection
        IpdXrealP = rxPdR .* abs(xHetP).^2;
        IpdXrealN = rxPdR .* abs(xHetN).^2;
        IpdYrealP = rxPdR .* abs(yHetP).^2;
        IpdYrealN = rxPdR .* abs(yHetN).^2;
        if sysParam.addThermNoise
            pd_thmvar = 4 * BOLTZMAN * TEMPERATURE / PD_LOAD_RESISTANCE * samplingFs;
            IpdThermal1 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_thmvar,'linear','real');
            IpdThermal2 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_thmvar,'linear','real');
            IpdThermal3 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_thmvar,'linear','real');
            IpdThermal4 = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_thmvar,'linear','real');
        else
            IpdThermal1      = 0;
            IpdThermal2      = 0;
            IpdThermal3      = 0;
            IpdThermal4      = 0;
        end
        
        if sysParam.addShotNoise
            pd_shtvar1   = 2*ELECTRON*(calcrms(xHetP)^2)*samplingFs; % assuming responsitivity is 1
            pd_shtvar2   = 2*ELECTRON*(calcrms(xHetN)^2)*samplingFs;
            pd_shtvar3   = 2*ELECTRON*(calcrms(yHetP)^2)*samplingFs;
            pd_shtvar4   = 2*ELECTRON*(calcrms(yHetN)^2)*samplingFs;
            IpdShot1      = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar1,'linear','real');
            IpdShot2      = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar2,'linear','real');
            IpdShot3      = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar3,'linear','real');
            IpdShot4      = genWGN(size(IpdXrealP,1),size(IpdXrealP,2),pd_shtvar4,'linear','real');
        else
            IpdShot1      = 0;
            IpdShot2      = 0;
            IpdShot3      = 0;
            IpdShot4      = 0;
        end
        
        V1      = IpdXrealP + IpdDark + IpdThermal1 + IpdShot1;
        V2      = IpdXrealN + IpdDark + IpdThermal2 + IpdShot2;
        V3      = IpdYrealP + IpdDark + IpdThermal3 + IpdShot3;
        V4      = IpdYrealN + IpdDark + IpdThermal4 + IpdShot4;
        
        if sysParam.doBalanced % balanced detection
            pd_xi			= V1 - V2;
            pd_yi			= V3 - V4;
        else
            pd_xi			= V1;
            pd_yi			= V3;
        end
        % low-pass filtering
        pd_xi           = real(ifft(fft(pd_xi).* rxPulseShapeFilter.freqRespBessel));
        pd_yi           = real(ifft(fft(pd_yi).* rxPulseShapeFilter.freqRespBessel));
    otherwise
        error(EID,'detection mode not supported !!');
end

if ctrlParam.doPlot
    if exist('FIG_EYE','var')
        plotEyeDiagram(FIG_EYE,pd_xi, baudrate, samplingFs, 'electrical');
    else
        FIG_EYE = plotEyeDiagram([],pd_xi, baudrate, samplingFs, 'electrical');
    end
    title('Real part of signal after photodetector');
%     figure(FIG_RECEIVED); plot(pd_xi,pd_xq,'.'); grid on
end


% add electrical spectrum
% keyboard;



%% ADC
% ad_head         = round(samplesPerSym/2);
ad_head         = 1;
ad_sps_sim      = samplesPerSym;

% need a realistic sampling rate
switch DETECTION_MODE
    case 'HOM'
        adc1 = pd_xi(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
        adc2 = pd_xq(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
        adc3 = pd_yi(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
        adc4 = pd_yq(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
    case 'HET'
        adc1 = pd_xi(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
        adc2 = pd_yi(ad_head:ad_sps_sim/ADC_SAMPLING_RATE:end);
    otherwise
        error(EID,'detection mode not supported !!');
end

% also need a ENOB

% add timing jitter here

if ctrlParam.doPlot
    figure(FIG_AFTERADC); plot(adc1(1:2:end),adc2(1:2:end),'.',adc1(2:2:end),adc2(2:2:end),'g.'); grid on;
    title('Signal after ADC');
end

% keyboard;

% dsp
if DSP_MODE == 0
    
    % go offline
	adc_out_len = length(adc1);
    
    % push only the center frame to the memory
	memory1 = [memory1; adc1(adc_out_len/FRAME_WINDOW+1 : end-adc_out_len/FRAME_WINDOW, 1)];
	memory2 = [memory2; adc2(adc_out_len/FRAME_WINDOW+1 : end-adc_out_len/FRAME_WINDOW, 1)];
	memory3 = [memory3; adc3(adc_out_len/FRAME_WINDOW+1 : end-adc_out_len/FRAME_WINDOW, 1)];
	memory4 = [memory4; adc4(adc_out_len/FRAME_WINDOW+1 : end-adc_out_len/FRAME_WINDOW, 1)];
	
	dspMemCount = dspMemCount + 1;
	
	if dspMemCount == DSO_MEMORY_LENGTH
        fprintf('In total %d frames are captured in the memory\n', dspMemCount);
		break;
    end
    
elseif DSP_MODE == 1
    
    % go real-time
    
	dspout1		= adc1 + 1j*adc2;
	dspout2		= adc3 + 1j*adc4;
	dspout1		= dspout1(1:2:end);
	dspout2		= dspout2(1:2:end);
end

%% Decision

if DSP_MODE == 1
	
	de_x = normalizeQam(dspout1,ALPHABET_SIZE);
	de_y = normalizeQam(dspout2,ALPHABET_SIZE);
	
	bit1 = slicerGrayQam(de_x,ALPHABET_SIZE);
	bit2 = slicerGrayQam(de_y,ALPHABET_SIZE);
	
	biterr = nnz(bit1(:,1:symbolsPerFrame)-refBitsX)+nnz(bit2(:,1:symbolsPerFrame)-refBitsY);
	
	fprintf('run # %d\t error count # %d\n',RUN,biterr);
end

end %% End of main loop

% this is the end of main loop <--------------
% this is the end of main loop <--------------
% this is the end of main loop <--------------

%% Offline DSP
% [dspout1,dspout2] = dspMain_built_150927(memory1,memory2,memory3,memory4,dspParam);
[dspout1,dspout2] = dspMain_built_151229(memory1,memory2,memory3,memory4,dspParam);

h1 = figure(34); plot(dspout1,'.'); grid on;
h2 = figure(35); plot(dspout2,'.'); grid on;
mngFigureWindow(h1,h2);



%% Decision
%
% de_x = normalizeQam(dspout1,ALPHABET_SIZE);
% de_y = normalizeQam(dspout2,ALPHABET_SIZE);
% 
% bit1 = slicerGrayQam(de_x,ALPHABET_SIZE);
% bit2 = slicerGrayQam(de_y,ALPHABET_SIZE);
% 
% biterr = nnz(bit1-bitRefx_offline)+nnz(bit2-bitRefy_offline);

return

