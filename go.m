function vM = go(vSet, LOG, VERBOSE, MODE)
% Fiber-optics with DSP
%
% datestr(datenum(now))
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

if nargin < 4
    MODE = 0; % 0 back-to-back
    % 1
    % 2
    % 3
end

if nargin < 3
    VERBOSE = 1;
end

if nargin < 2
    LOG = 0;
end

% Log file
if LOG
    logFile = fopen('log','a');
    fprintf(logFile, '-------------------------------------');
    fprintf(logFile, '\n');
    fprintf(logFile, 'Time now is %s\n', datestr(now));
    if exist('vSet', 'var')
        fprintf(logFile, 'Caller is %s\n', vSet.caller);
    end
else
    logFile = 1;
end

switch MODE
    case 0
        ctrlParam.doPilot           = 0;
        ctrlParam.doNFC             = 0;
        ctrlParam.doRZ              = 0;
        ctrlParam.doNyquist         = 1;
        ctrlParam.doRndPMD          = 0;
        ctrlParam.doMzmComp         = 1;
        % ctrlParam.doCoherent        = 1;
        ctrlParam.doBalanced        = 1;
        ctrlParam.doDSP             = 1;
        ctrlParam.doPlot            = 0;
        % Controller
        ctrlParam.addASE            = 1;
        ctrlParam.addCD             = 0;
        ctrlParam.addPMD            = 0;
        ctrlParam.addLaserRIN       = 0;
        ctrlParam.addLaserPN        = 0;
        ctrlParam.addPolarRot       = 0;
        ctrlParam.addThermNoise     = 0;
        ctrlParam.addShotNoise      = 0;
        ctrlParam.addRIN            = 0;
        ctrlParam.addAdcClockJitter = 0;
        ctrlParam.addFreqOffset     = 0;
    case 1
    case 2
    case 3
    otherwise
        keyboard;
end

% Constants
LIGHT_SPEED             = 299792458;
BOLTZMAN                = 1.381e-23;
ELECTRON                = 1.602e-19;
EID                     = 'goErr';
FRAME_WINDOW            = 3;
TEMPERATURE             = 300;
PD_LOAD_RESISTANCE      = 5000;             % TIA in PD
NOISE_REFERENCE_BAND    = 12.5e9;           % for OSNR definition
DETECTION_MODE          = 'HOM';            % HOM or HET
CENTER_FREQUENCY        = 193.41e12;
CENTER_WAVELENGTH       = LIGHT_SPEED / CENTER_FREQUENCY;
FIG_TXPN                = 1;
FIG_TXLASER             = 2;
FIG_DRVEYE              = 3;
FIG_RECEIVED            = 4;
FIG_AFTERADC            = 5;

if nargin < 1
    MAX_RUN_NUMBER          = 1000;
    HYBRID_90_PHASESHIFT    = 90;           % degree
    ADC_SAMPLING_RATE       = 2;            % samples per symbol
    DSP_MODE                = 0;            % 0 - offline; 
    DSO_MEMORY_LENGTH       = 200;          % number of frames
    LASER_LINEWIDTH         = 500e3;
    OSNR                    = 140;           % SNR of one symbol
    baudrate                = 30e9;
    bitpersym               = 2;
    modFormat               = 'QPSK';
    freqOffset              = 1.0e9;
    psFiltType              = 'Nyquist';    % Nyquist Bessel Gaussian
else
    MAX_RUN_NUMBER          = vSet.nFrm;
    HYBRID_90_PHASESHIFT    = vSet.Hyd90;
    ADC_SAMPLING_RATE       = vSet.ADCfs;
    DSP_MODE                = vSet.DSPmode;
    DSO_MEMORY_LENGTH       = vSet.DSPmemLen; 
    LASER_LINEWIDTH         = vSet.linewidth;
    OSNR                    = vSet.osnr;
    baudrate                = vSet.baudrate;
    bitpersym               = vSet.bitpersym;
    modFormat               = vSet.modFormat;
    freqOffset              = vSet.freqOffset;
    psFiltType              = vSet.psFiltType;
end

if VERBOSE
    fprintf(logFile, '* modulation format is %s \n', modFormat);
    fprintf(logFile, '* pulse-shaping filter type is %s \n', psFiltType);
    fprintf(logFile, '* detection mode is %s \n', DETECTION_MODE);
end

% Global parameters
samplingFs          = 8 * baudrate;
timewindow          = 512 / baudrate;
symbolsPerFrame     = timewindow * baudrate;
samplePerSymbol     = samplingFs / baudrate;
samplesPerFrame     = symbolsPerFrame * samplePerSymbol;
vctFreqPerFrm       = getFFTGrid(samplesPerFrame, samplingFs);
numSamples          = FRAME_WINDOW * samplesPerFrame;
timeVector          = (0 : numSamples-1)' / samplingFs;
freqVector          = getFFTGrid(numSamples, samplingFs);
ALPHABET_SIZE       = 2 ^ bitpersym;
vM.StartTime        = datestr(now);

% initialize absolute time
timeVectorAbs       = timeVector;

% Components
txLaser.centerFreq = CENTER_FREQUENCY;
txLaser.centerLambda = LIGHT_SPEED / txLaser.centerFreq;
txLaser.linewidth = LASER_LINEWIDTH;
txLaser.phaseNoiseVar = 2 * pi * txLaser.linewidth / samplingFs;
txLaser.azimuth = 0;
txLaser.ellipticity = 0;
txLaser.power = 1e-3;
txLaser.RIN = -130;                     % dBc/Hz

modulator.Vpi = 3;                      % [V]
modulator.extRatio = 100;               % [dB]
modulator.efficiency = 0.15;

transmtter.lpfOrder = 4;
transmtter.bandwidth = 0.75 * baudrate;
transmtter.pilotX = [];
transmtter.pilotY = [];

fiber.n2 = 2.6e-20;
fiber.coreArea = 80e-12;                % [m^2]
fiber.initPolAngle = 0;                 % degree
fiber.rotSpeed = 30e3;                  % rad/s
fiber.lossFast = 0.2;
fiber.lossSlow = 0.2;                   % [dB/km]
fiber.spanLength = 80e3;
fiber.spanNum = 2;
fiber.stepLength = 1e3;
fiber.corrLength = 100;
fiber.dispParamD = 17e-6;               % [s/m]
fiber.dispParamS = 0.08e3;              % [s/m^2]
fiber.doFullPMD = 0;
fiber.pmdParam = 0.5e-12 / 31.623;
fiber.noiseSeed = 0;                    % random number seed

% accumulated dispersion ps/nm
fiber.DL = fiber.dispParamD * fiber.spanLength * fiber.spanNum;
% accumulated dispersion slop
fiber.SL = fiber.dispParamS * fiber.spanLength * fiber.spanNum;

optFilterOrder = 4;
optFilterBw = 40e9;

receiver.centerFreq = CENTER_FREQUENCY;
receiver.lpfOrder = 5;
receiver.bandwidth = 0.75 * baudrate;

rxLaser.centerLambda = LIGHT_SPEED / receiver.centerFreq;
rxLaser.noiseSeed = 0;
rxLaser.linewidth = LASER_LINEWIDTH;
rxLaser.phaseNoiseVar = 2 * pi * rxLaser.linewidth / samplingFs;
rxLaser.azimuth = 45;
rxLaser.ellipticity = 0;
rxLaser.power = 10e-3;
rxLaser.psdRIN = -135;                  % dBc/Hz, one-sided

% power split ratio of PBC
receiver.pbcPowerSplitRatio   = 0.5;
receiver.pbcPhaseRetard    = 0;
receiver.detectorResponsivity = 1.0;
receiver.detectorDarkCurrent = 0.0E-9;
receiver.CMRR   = 30;

% electrical noises are defined as one-sided PSD in current
receiver.psdThermal = 4 * BOLTZMAN * TEMPERATURE / PD_LOAD_RESISTANCE * (0.5 * samplingFs);
receiver.psdShot    = 2 * ELECTRON * receiver.detectorResponsivity * (rxLaser.power / 2) * (0.5 * samplingFs);
receiver.psdRIN     = (receiver.detectorResponsivity).^2 * (10.^(-receiver.CMRR/20)).^2 * idbw(rxLaser.psdRIN) * (rxLaser.power).^2 * (0.5 * samplingFs) / 4;


%% Preparing filter responses
txPulseShapeFilter.RollOffFactor = 0.35;
txPulseShapeFilter.freqRespRC = calcFilterFreqResp(numSamples, samplingFs, txPulseShapeFilter.RollOffFactor, baudrate, 'rc');
txPulseShapeFilter.freqRespRRC = calcFilterFreqResp(numSamples, samplingFs, txPulseShapeFilter.RollOffFactor, baudrate, 'rrc');

txPulseShapeFilter.freqRespBessel = calcFilterFreqResp(numSamples, samplingFs, transmtter.lpfOrder, transmtter.bandwidth, 'bessel');
txPulseShapeFilter.freqRespNyquist = calcFilterFreqResp(numSamples, samplingFs, 0.1, baudrate, 'nyquist');
txPulseShapeFilter.freqRespGaussian = calcFilterFreqResp(numSamples, samplingFs, transmtter.lpfOrder, transmtter.bandwidth, 'gaussian');


rxPulseShapeFilter.RollOffFactor = 0.35;
rxPulseShapeFilter.freqRespRC = calcFilterFreqResp(numSamples, samplingFs, rxPulseShapeFilter.RollOffFactor, baudrate, 'rc');
rxPulseShapeFilter.freqRespRRC = calcFilterFreqResp(numSamples, samplingFs, rxPulseShapeFilter.RollOffFactor, baudrate, 'rrc');

rxPulseShapeFilter.freqRespBessel = calcFilterFreqResp(numSamples, samplingFs, receiver.lpfOrder, receiver.bandwidth, 'bessel');
rxPulseShapeFilter.freqRespNyquist = calcFilterFreqResp(numSamples, samplingFs, 0.1, baudrate, 'nyquist');
rxPulseShapeFilter.freqRespGaussian = calcFilterFreqResp(numSamples, samplingFs, receiver.lpfOrder, receiver.bandwidth, 'gaussian');

optGaussianFilter = calcOptGaussFlt(numSamples, samplingFs, 0, optFilterOrder, optFilterBw);


% Initialize buffer
buffer.txBitsX = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.txBitsY = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.rxBitsX = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.rxBitsY = randi([0 1], bitpersym, FRAME_WINDOW*symbolsPerFrame);
buffer.txPhaseNoise = genLaserPhaseNoise(numSamples, txLaser.phaseNoiseVar, 0);
buffer.rxPhaseNoise = genLaserPhaseNoise(numSamples, rxLaser.phaseNoiseVar, 0);
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

% Clear memory
memory1 = [];
memory2 = [];
memory3 = [];
memory4 = [];

dspMemCount = 0;

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
dspParam.doCPE				= 0;
dspParam.doLmsAfterCPE      = 0;
dspParam.doFOE				= 0;
dspParam.doFEC				= 0;
dspParam.doDownSampling     = 1;

dspParam.DL = fiber.DL;
dspParam.SL = fiber.SL;
dspParam.lambda = CENTER_WAVELENGTH;
dspParam.lambda0 = CENTER_WAVELENGTH;

dspParam.lpfBW = 0;

dspParam.mimoStep = 1e-3;
dspParam.mimoTapNum = 7;                %7 or 13
dspParam.mimoErrId = 0;
dspParam.mimoIteration = 10;

dspParam.lmsGainAfterCPE = 1e-3;
dspParam.lmsTapsAfterCPE = 125;
dspParam.lmsIterAfterCPE = 4;

dspParam.cpeBlockSize = 16;
dspParam.cpePLLmu = [1e-3, 1e-6];
dspParam.doMLCPE = 2;
dspParam.cpeMLavgLeng = 32;
dspParam.cpeMlIter = 1;
dspParam.vvpeAvgMode = 1;               % 0 for blockwise 1 for sliding window
dspParam.cpeBPSnTestPhase = 10;
dspParam.cpeAlgSelect = 'VVPE';

refBitsOfflineX		= [];
refBitsOfflineY		= [];

% Start main loop
for RUN = 1 : MAX_RUN_NUMBER
    % Binary source for a new frame
    bitsX = randi([0 1],bitpersym,symbolsPerFrame);
    bitsY = randi([0 1],bitpersym,symbolsPerFrame);
    
    % fifo
    [buffer.txBitsX] = fifoBuffer(buffer.txBitsX,bitsX);
    [buffer.txBitsY] = fifoBuffer(buffer.txBitsY,bitsY);
    
    % offline mode; take the center frame as the reference
    if DSP_MODE == 0
        refBitsOfflineX = [refBitsOfflineX buffer.txBitsX(:,(1 : symbolsPerFrame) + symbolsPerFrame)];
        refBitsOfflineY = [refBitsOfflineY buffer.txBitsY(:,(1 : symbolsPerFrame) + symbolsPerFrame)];
    elseif DSP_MODE == 1
        refBitsX = buffer.txBitsX(:,1 : symbolsPerFrame);
        refBitsY = buffer.txBitsY(:,1 : symbolsPerFrame);
    end
    
    % Bit mapping
    % convert bits to syms
    txBaudX = symbolizerGrayQam(buffer.txBitsX);
    txBaudY = symbolizerGrayQam(buffer.txBitsY);
    
    % polarizationAnalyzer(txBaudX.',txBaudY.','o-');
    
    % TX laser 
    % generate new phase noise for the new frame
    tmpPN = genLaserPhaseNoise(samplesPerFrame, txLaser.phaseNoiseVar, buffer.txPhaseNoise(end));
    
    % fifo
    buffer.txPhaseNoise = fifoBuffer(buffer.txPhaseNoise, tmpPN);
    
    % buffer laser RIN is really defined as one-sided power variance
    
    % generate new rin for new frame
    tmpRIN = genWGN(1, samplesPerFrame, idbw(txLaser.RIN) * (txLaser.power.^2) * (0.5 * samplingFs), 'linear', 'real');
    
    % fifo
    buffer.txLaserRIN = fifoBuffer(buffer.txLaserRIN, tmpRIN);
    
    if ctrlParam.doPlot
        pltIndex = 1:500;
        figure(FIG_TXPN); plot(tmpPN(pltIndex)); grid on; title('Tx Phase Noise');
    end
    
    % laser on
    if ctrlParam.addLaserRIN
        if ctrlParam.addLaserPN
            txLaser.wfm = sqrt(txLaser.power + buffer.txLaserRIN) .* exp(1j * buffer.txPhaseNoise);
        else
            txLaser.wfm = sqrt(txLaser.power + buffer.txLaserRIN);
        end
    else
        if ctrlParam.addLaserPN
            txLaser.wfm = sqrt(txLaser.power) .*  exp(1j * buffer.txPhaseNoise);
        else
            txLaser.wfm = sqrt(txLaser.power) * ones(1,numSamples);
        end
    end
    
    % if ctrlParam.doPlot
    %     figure(FIG_TXLASER); plot(abs(txLaser(pltIndex).^2)); grid on; title('Tx laser power waveform');
    % end
    
    % Driver
    % normalize
    txBaudRealX = real(txBaudX) / (sqrt(ALPHABET_SIZE)-1);
    txBaudImagX = imag(txBaudX) / (sqrt(ALPHABET_SIZE)-1);
    txBaudRealY = real(txBaudY) / (sqrt(ALPHABET_SIZE)-1);
    txBaudImagY = imag(txBaudY) / (sqrt(ALPHABET_SIZE)-1);
    
    % Pulse shaping
    % DAC - simple oversampling by inserting zeros
    txDrvIxUps = upSampInsertZeros(txBaudRealX, samplePerSymbol);
    txDrvQxUps = upSampInsertZeros(txBaudImagX, samplePerSymbol);
    txDrvIyUps = upSampInsertZeros(txBaudRealY, samplePerSymbol);
    txDrvQyUps = upSampInsertZeros(txBaudImagY, samplePerSymbol);
    
    % filtering
    txDrvIxWfm = real(ifft(fft(txDrvIxUps) .* txPulseShapeFilter.freqRespRRC));
    txDrvQxwfm = real(ifft(fft(txDrvQxUps) .* txPulseShapeFilter.freqRespRRC));
    txDrvIyWfm = real(ifft(fft(txDrvIyUps) .* txPulseShapeFilter.freqRespRRC));
    txDrvQyWfm = real(ifft(fft(txDrvQyUps) .* txPulseShapeFilter.freqRespRRC));
    
%     plotEyeDiagram(txDrvIxWfm, 2*samplingFs/baudrate, 'e');

    % pre-distortion
    if ctrlParam.doMzmComp
        txDrvIxWfm = asin(txDrvIxWfm) * modulator.Vpi /pi;
        txDrvQxwfm = asin(txDrvQxwfm) * modulator.Vpi /pi;
        txDrvIyWfm = asin(txDrvIyWfm) * modulator.Vpi /pi;
        txDrvQyWfm = asin(txDrvQyWfm) * modulator.Vpi /pi;
    else
        txDrvIxWfm = (txDrvIxWfm * modulator.efficiency) * modulator.Vpi /pi;
        txDrvQxwfm = (txDrvQxwfm * modulator.efficiency) * modulator.Vpi /pi;
        txDrvIyWfm = (txDrvIyWfm * modulator.efficiency) * modulator.Vpi /pi;
        txDrvQyWfm = (txDrvQyWfm * modulator.efficiency) * modulator.Vpi /pi;
    end
    
    % MZM non-linear pre-comp
    % *add MZM nonlinearity pre-comp. here...*
    
    % MZM
    V1 = - modulator.Vpi / 2;
    V2 = - modulator.Vpi / 2;
    V3 = + modulator.Vpi / 2;
    % the nonlinear transfer function of MZM will degradate the matched
    % filtering
    txOptSigX = oeModIqNested(txLaser.wfm, txDrvIxWfm, txDrvQxwfm, modulator.extRatio, modulator.Vpi, V1, V2, V3);
    txOptSigY = oeModIqNested(txLaser.wfm, txDrvIyWfm, txDrvQyWfm, modulator.extRatio, modulator.Vpi, V1, V2, V3);
    
%     plotEyeDiagram(txOptSigX, 2*samplingFs/baudrate, 'o');
    
    % OSNR
    % here goes an OSNR emulator
    
    % Fiber channel
    if ~ ctrlParam.doRndPMD % if random birefrigence is switched off, use simple model
        % convert osnr to snr per sample (oversampling)
        SNR = osnr2snr(OSNR, baudrate, samplePerSymbol, 'complex');
        % snr per symbol is much higher than snr oversampling
        vM.SNR = SNR + dbw(samplePerSymbol);
        vM.OSNR = OSNR;
        
        % calc noise power, pulse-shaping is power perserving
        sigPowXdb   = dbw(calcrms(txOptSigX).^2);
        sigPowYdb   = dbw(calcrms(txOptSigY).^2);
        noisePowerx = idbw(sigPowXdb - SNR);
        noisePowery = idbw(sigPowYdb - SNR);
        
        % generate new noise for new frame
        channelNoiseX = genWGN(1, samplesPerFrame, noisePowerx, 'linear', 'complex');
        channelNoiseY = genWGN(1, samplesPerFrame, noisePowery, 'linear', 'complex');
        
        % fifo
        buffer.channelNoiseX = fifoBuffer(buffer.channelNoiseX, channelNoiseX);
        buffer.channelNoiseY = fifoBuffer(buffer.channelNoiseY, channelNoiseY);
        
        % add noise
        if ctrlParam.addASE
            txOptSigX = txOptSigX + buffer.channelNoiseX(:);
            txOptSigY = txOptSigY + buffer.channelNoiseY(:);
        end
        
        % fiber
        lk_rotjump = mod(fiber.rotSpeed * timeVectorAbs, 2*pi);
        lk_theta = fiber.initPolAngle + lk_rotjump;
        %     fiber.initPolAngle = lk_theta(numSamples);
        
        if ctrlParam.addPolarRot % rotate device only if polarization rotation is on
            tmpOptSigX = txOptSigX .* cos(lk_theta) - txOptSigY .* sin(lk_theta);
            tmpOptSigY = txOptSigX .* sin(lk_theta) + txOptSigY .* cos(lk_theta);
        else
            tmpOptSigX = txOptSigX;
            tmpOptSigY = txOptSigY;
        end
        
        if ctrlParam.addCD
            [HCD] = calcDispResponse(numSamples, samplingFs, centerWave, centerWave, fiber.DL, fiber.SL);
            tmpOptSigX = ifft(fft(tmpOptSigX) .* HCD);
            tmpOptSigY = ifft(fft(tmpOptSigY) .* HCD);
        end
        
        if ctrlParam.addPMD
            tmpOptSigX = ifft(fft(tmpOptSigX) .* exp(-1i * 2 * pi * freqVector * (+ lk_DGD / 2)) .* exp(-1i * 2 * pi * uv_fc * (+ lk_DGD / 2)));
            tmpOptSigY = ifft(fft(tmpOptSigY) .* exp(-1i * 2 * pi * freqVector * (- lk_DGD / 2)) .* exp(-1i * 2 * pi * uv_fc * (- lk_DGD / 2)));
        end
        
        if ctrlParam.addPolarRot % rotate back only if polarization rotation is on
            tmpOptSigX = tmpOptSigX.*cos(-lk_theta) - tmpOptSigY.*sin(-lk_theta);
            tmpOptSigY = tmpOptSigX.*sin(-lk_theta) + tmpOptSigY.*cos(-lk_theta);
        end
        
        if ctrlParam.addPolarRot % add endless rotation
            tmpOptSigX = tmpOptSigX .* cos(lk_theta) - tmpOptSigY .* sin(lk_theta);
            tmpOptSigY = tmpOptSigX .* sin(lk_theta) + tmpOptSigY .* cos(lk_theta);
            % update absolute time
            timeVectorAbs = ((0 : numSamples-1)' + samplesPerFrame) / samplingFs;
            lk_rotjump = mod(fiber.rotSpeed * timeVectorAbs, 2*pi);
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
        spectrumAnalyzer(txOptSigX, freqVector);
        spectrumAnalyzer(rxOptSigX, freqVector);
    end
    
    
    % RX laser
    % polarization control
    loJonesVector = calcJonesVector(rxLaser.azimuth, rxLaser.ellipticity);
    
    % generate new phase noise for new frame
    tmpPN = genLaserPhaseNoise(samplesPerFrame, rxLaser.phaseNoiseVar, buffer.rxPhaseNoise(end));
    
    % fifo
    buffer.rxPhaseNoise = fifoBuffer(buffer.rxPhaseNoise, tmpPN);
    
    % buffer laser rin
    tmpRIN = genWGN(1, samplesPerFrame, idbw(rxLaser.psdRIN) * (rxLaser.power.^2) * (samplingFs / 2), 'linear', 'real');
    buffer.rxLaserRIN = fifoBuffer(buffer.rxLaserRIN, tmpRIN);
    
    % laser on
    if ctrlParam.addLaserRIN
        if ctrlParam.addLaserPN
            rxLaser.wfm = sqrt(rxLaser.power + buffer.rxLaserRIN) .* exp(1i * buffer.rxPhaseNoise);
        else
            rxLaser.wfm = sqrt(rxLaser.power + buffer.rxLaserRIN);
        end
    else
        if ctrlParam.addLaserPN
            rxLaser.wfm = sqrt(rxLaser.power) .* exp(1i * buffer.rxPhaseNoise);
        else
            rxLaser.wfm = sqrt(rxLaser.power) * ones(1, length(buffer.rxLaserRIN));
        end
    end
    
    % control the polarization
    rxLaser.wfm = loJonesVector * sqrt(abs(rxLaser.wfm).^2) .* [sign(rxLaser.wfm); sign(rxLaser.wfm)];
    
    % PBS / PBC
    loLaserPx = sqrt(receiver.pbcPowerSplitRatio) * exp(-1i * receiver.pbcPhaseRetard) * rxLaser.wfm(1,:).';
    loLaserPy = sqrt(1 - receiver.pbcPowerSplitRatio) * rxLaser.wfm(2,:).';
    
    
    % Hybrid
    hybrid90 = deg2rad(HYBRID_90_PHASESHIFT);
    
    switch DETECTION_MODE
        case 'HOM'
            if ctrlParam.addFreqOffset
                xRealP = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .*  loLaserPx;
                xRealN = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .* -loLaserPx;
                xImagP = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .*  exp(1i * hybrid90) .* loLaserPx;
                xImagN = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .* -exp(1i * hybrid90) .* loLaserPx;
                yRealP = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .*  loLaserPy;
                yRealN = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .* -loLaserPy;
                yImagP = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .*  exp(1i * hybrid90) .* loLaserPy;
                yImagN = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .* -exp(1i * hybrid90) .* loLaserPy;
            else
                xRealP = rxOptSigX + loLaserPx;
                xRealN = rxOptSigX - loLaserPx;
                xImagP = rxOptSigX + exp(1i * hybrid90) .* loLaserPx;
                xImagN = rxOptSigX - exp(1i * hybrid90) .* loLaserPx;
                yRealP = rxOptSigY + loLaserPy;
                yRealN = rxOptSigY - loLaserPy;
                yImagP = rxOptSigY + exp(1i * hybrid90) .* loLaserPy;
                yImagN = rxOptSigY - exp(1i * hybrid90) .* loLaserPy;
            end
        case 'HET'
            if ctrlParam.addFreqOffset
                xHetP = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .*  loLaserPx;
                xHetN = rxOptSigX + exp(1i * 2*pi * freqOffset * timeVector) .* -loLaserPx;
                yHetP = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .*  loLaserPy;
                yHetN = rxOptSigY + exp(1i * 2*pi * freqOffset * timeVector) .* -loLaserPy;
            else
                xHetP = rxOptSigX + loLaserPx;
                xHetN = rxOptSigX - loLaserPx;
                yHetP = rxOptSigY + loLaserPy;
                yHetN = rxOptSigY - loLaserPy;
            end
        otherwise
            warning('detection mode not supported'); keyboard;
    end
    
    
    % Photo detector
    % dark current
    IpdDark = receiver.detectorDarkCurrent;
    
    switch DETECTION_MODE
        case 'HOM'
            % square law detection
            V1 = receiver.detectorResponsivity .* abs(xRealP).^2;
            V2 = receiver.detectorResponsivity .* abs(xRealN).^2;
            V3 = receiver.detectorResponsivity .* abs(xImagP).^2;
            V4 = receiver.detectorResponsivity .* abs(xImagN).^2;
            V5 = receiver.detectorResponsivity .* abs(yRealP).^2;
            V6 = receiver.detectorResponsivity .* abs(yRealN).^2;
            V7 = receiver.detectorResponsivity .* abs(yImagP).^2;
            V8 = receiver.detectorResponsivity .* abs(yImagN).^2;
            
            if ctrlParam.addThermNoise
                IpdThermal1 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
                IpdThermal2 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
                IpdThermal3 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
                IpdThermal4 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
            else
                IpdThermal1 = 0;
                IpdThermal2 = 0;
                IpdThermal3 = 0;
                IpdThermal4 = 0;
            end
            
            if ctrlParam.addShotNoise
                IpdShot1 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
                IpdShot2 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
                IpdShot3 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
                IpdShot4 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
            else
                IpdShot1 = 0;
                IpdShot2 = 0;
                IpdShot3 = 0;
                IpdShot4 = 0;
            end
            
            if ctrlParam.addRIN
                IpdRIN1 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
                IpdRIN2 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
                IpdRIN3 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
                IpdRIN4 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
            else
                IpdRIN1 = 0;
                IpdRIN2 = 0;
                IpdRIN3 = 0;
                IpdRIN4 = 0;
            end
            
            if ctrlParam.doBalanced 
                pd_xi = V1 - V2 + IpdThermal1 + IpdShot1 + IpdDark + IpdRIN1;
                pd_xq = V3 - V4 + IpdThermal2 + IpdShot2 + IpdDark + IpdRIN2;
                pd_yi = V5 - V6 + IpdThermal3 + IpdShot3 + IpdDark + IpdRIN3;
                pd_yq = V7 - V8 + IpdThermal4 + IpdShot4 + IpdDark + IpdRIN4;
            else
                pd_xi = V1 + IpdThermal1 + IpdShot1 + IpdDark;
                pd_xq = V3 + IpdThermal2 + IpdShot2 + IpdDark;
                pd_yi = V5 + IpdThermal3 + IpdShot3 + IpdDark;
                pd_yq = V7 + IpdThermal4 + IpdShot4 + IpdDark;
            end
            
            % match filtering
            pd_xi = real(ifft(fft(pd_xi) .* rxPulseShapeFilter.freqRespRRC));
            pd_xq = real(ifft(fft(pd_xq) .* rxPulseShapeFilter.freqRespRRC));
            pd_yi = real(ifft(fft(pd_yi) .* rxPulseShapeFilter.freqRespRRC));
            pd_yq = real(ifft(fft(pd_yq) .* rxPulseShapeFilter.freqRespRRC));
            
        case 'HET'
            % square law detection
            V1 = receiver.detectorResponsivity .* abs(xHetP).^2;
            V2 = receiver.detectorResponsivity .* abs(xHetN).^2;
            V5 = receiver.detectorResponsivity .* abs(yHetP).^2;
            V6 = receiver.detectorResponsivity .* abs(yHetN).^2;
            if ctrlParam.addThermNoise
                IpdThermal1 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
                IpdThermal2 = genWGN(size(V1,1), size(V1,2), receiver.psdThermal, 'linear', 'real');
            else
                IpdThermal1 = 0;
                IpdThermal2 = 0;
            end
            
            if ctrlParam.addShotNoise
                IpdShot1 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
                IpdShot2 = genWGN(size(V1,1), size(V1,2), receiver.psdShot, 'linear', 'real');
            else
                IpdShot1      = 0;
                IpdShot2      = 0;
            end
            
            if ctrlParam.addRIN
                IpdRIN1 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
                IpdRIN2 = genWGN(size(V1,1), size(V1,2), receiver.psdRIN, 'linear', 'real');
            else
                IpdRIN1 = 0;
                IpdRIN2 = 0;
            end

            if ctrlParam.doBalanced % balanced detection
                pd_xi = V1 - V2 + IpdDark + IpdThermal1 + IpdShot1 + IpdRIN1;
                pd_yi = V3 - V4 + IpdDark + IpdThermal2 + IpdShot2 + IpdRIN2;
            else
                pd_xi = V1 + IpdDark + IpdThermal1 + IpdShot1;
                pd_yi = V3 + IpdDark + IpdThermal2 + IpdShot2;
            end
            
            % low-pass filtering
            pd_xi = real(ifft(fft(pd_xi) .* rxPulseShapeFilter.freqRespBessel));
            pd_yi = real(ifft(fft(pd_yi) .* rxPulseShapeFilter.freqRespBessel));
            
        otherwise
            warning('detection mode not supported'); keyboard;
    end
    
    if ctrlParam.doPlot
        plotEyeDiagram(pd_xi, samplingFs / baudrate, 'e');
        title('Real part of signal after photodetector');
    end
    
    % add electrical spectrum
    % keyboard;
    
    % ADC
    % ad_head         = round(samplePerSymbol/2);
    ad_head         = 1;
    ad_sps_sim      = samplePerSymbol;
    
    % need a realistic sampling rate
    switch DETECTION_MODE
        case 'HOM'
            adc1 = pd_xi(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
            adc2 = pd_xq(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
            adc3 = pd_yi(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
            adc4 = pd_yq(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
        case 'HET'
            adc1 = pd_xi(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
            adc2 = pd_yi(ad_head : ad_sps_sim/ADC_SAMPLING_RATE : end);
        otherwise
            warning('detection mode not supported'); keyboard;
    end
    
    % also need a ENOB
    
    % add timing jitter here
    
    % dsp
    if DSP_MODE == 0
        %
        % go offline
        %
        adc_out_len = length(adc1);
        
        % push only the center frame to the memory
        memory1 = [memory1; adc1(adc_out_len / FRAME_WINDOW + 1 : end - adc_out_len / FRAME_WINDOW, 1)];
        memory2 = [memory2; adc2(adc_out_len / FRAME_WINDOW + 1 : end - adc_out_len / FRAME_WINDOW, 1)];
        memory3 = [memory3; adc3(adc_out_len / FRAME_WINDOW + 1 : end - adc_out_len / FRAME_WINDOW, 1)];
        memory4 = [memory4; adc4(adc_out_len / FRAME_WINDOW + 1 : end - adc_out_len / FRAME_WINDOW, 1)];
        
        dspMemCount = dspMemCount + 1;
        
        if dspMemCount == DSO_MEMORY_LENGTH
            fprintf('In total %d frames are captured in the memory\n', dspMemCount);
            break;
        end
        
    elseif DSP_MODE == 1
        %
        % go real-time
        %
        dspout1 = adc1 + 1i * adc2;
        dspout2 = adc3 + 1i * adc4;
        dspout1 = dspout1(1 : 2 : end);
        dspout2 = dspout2(1 : 2 : end);
    end
    
    % Decision
    if DSP_MODE == 1
        
        de_x = normalizeQam(dspout1, ALPHABET_SIZE);
        de_y = normalizeQam(dspout2, ALPHABET_SIZE);
        
        bit1 = slicerGrayQam(de_x, ALPHABET_SIZE);
        bit2 = slicerGrayQam(de_y, ALPHABET_SIZE);
        
        biterr = nnz(bit1(:, 1 : symbolsPerFrame) - refBitsX) + nnz(bit2(:, 1 : symbolsPerFrame) - refBitsY);
        
        fprintf('run # %d\t error count # %d\n', RUN, biterr);
    end
end %% End of main loop

% this is the end of main loop <--------------
% this is the end of main loop <--------------
% this is the end of main loop <--------------


% Offline DSP
if VERBOSE
    fprintf('Starting DSP...\n');
end
[dspout1, dspout2] = dspMain_built_151229(memory1, memory2, memory3, memory4, dspParam);

h1 = figure(34); plot(dspout1,'.'); grid on;
h2 = figure(35); plot(dspout2,'.'); grid on;
mngFigureWindow(h1, h2);

bit1 = slicerGrayQam(dspout1, ALPHABET_SIZE);
bit2 = slicerGrayQam(dspout2, ALPHABET_SIZE);

berx = sum(sum(abs(bit1 - refBitsOfflineX))) / numel(bit1);
bery = sum(sum(abs(bit2 - refBitsOfflineY))) / numel(bit2);

if VERBOSE
    fprintf('BER x = %.2e \n', berx);
    fprintf('BER y = %.2e \n', bery);
end

vM.BER = 0.5 * (berx + bery);
vM.BER_Theo = snr2ber(vM.SNR, bitpersym, 'dB');
if vM.BER < vM.BER_Theo
    warning('SANITY CHECK: Break BER limit, %.2E (%.2E in theory)', vM.BER, vM.BER_Theo);
end
return

