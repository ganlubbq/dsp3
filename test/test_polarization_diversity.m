%% Single Run Test Script for Maximum ratio combining (MRC)
%
% MRC is the optimal solution of utilizing the receiver diversity in terms
% of SNR;
%
% Since the data is transmitted in single polarization, the receiver
% firstly estimate the rotation matrix M of the polar-diversity coherent
% receiver including the power split ratio and phase difference and then
% rotate the signal back into single poalrization;
% 
% The improvment on average SNR is linearly proportional to # of receivers;
%
% Simple explanation is the power of MRC signal is doubled while the
% variance of noise remains the same;
%
% MIMO structured equalization has comparable (better?) performance than
% the MRC;
%
% Maximum Ratio Combining: multiply each received data with conjugate of
% its channel parameter
%
% Equalize Gain Combining: do not rotate, just normalize and add;
%
% Note on Receive diversity by Prof. RaviRaj Adve.
% http://www.comm.utoronto.ca/~rsadve/Notes/DiversityReceive.pdf

clear 
clc

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));

%% Global setting

%-----------constant
LIGHT_SPEED             = 299792458;
BOLTZMAN                = 1.381e-23;
ELECTRON                = 1.602e-19;
MAX_RUN_NUMBER          = 1000;
MAX_FIGURE_NUMBER       = 50;
FRAME_WINDOW            = 3;
TEMPERATURE             = 300;
PD_LOAD_RESISTANCE      = 5000; % TIA in PD
NOISE_REFERENCE_BAND    = 12.5e9;
HYBRID_90_PHASESHIFT    = 90; % degree
ADC_SAMPLING_RATE       = 2; % samples per symbol
DSP_MODE                = 0; % 0-offline; 1-real time
DSO_MEMORY_LENGTH       = 64; % number of frames

LASER_LINEWIDTH		= 0;
OSNR                = 15;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symlen = 2^12;
symrate = 30e9;

swing = [0 1];
bps = 2;

ALPHABET_SIZE = 2^bps;

% generate new stuff
bits = randi(swing,bps,symlen);


%% Bit mapping

% convert bits to syms
syms = graySymbolizer_mQAM(bits);

deltaLsr = LASER_LINEWIDTH/symrate;
lsrPNvar = 2*pi*deltaLsr;
lsrPN = genLaserPhaseNoise(symlen,lsrPNvar,0);
symLsr = syms .* exp(1j*lsrPN);

EsNo = OSNR - 10*log10(symrate/NOISE_REFERENCE_BAND);

Psig = dbw((calcrms(symLsr))^2);

Pns = Psig - EsNo;


%% PBS

% power split ratio
alpha = 0.55;

% phase retardation
delta = 0.1*pi;

H_PBS = [sqrt(alpha)*exp(1j*delta); sqrt(1-alpha)];

symOut = H_PBS * symLsr.';


%% Add noise

symTx(1,:) = symOut(1,:) + genWGN(length(symOut),Pns,'dbw','complex');
symTx(2,:) = symOut(2,:) + genWGN(length(symOut),Pns,'dbw','complex');
% scatterplot(symTx);


%% Verify the MRC

r = mean(symTx(1,:)./symTx(2,:));

if abs(r) < 1
    delta = angle(r);
    alpha = abs(r)^2 / (abs(r)^2 + 1);
else
    delta = -angle(1/r);
    alpha = 1 / (abs(1/r)^2 + 1);
end

M = [sqrt(alpha)*exp(-1j*delta), sqrt(1-alpha); ...
     -sqrt(1-alpha)*exp(-1j*delta), sqrt(alpha)]
 
% the first row of M is conjugate of H_PBS
 
SYMOUT_MRC = M * symTx;

mngScatterplot(symTx(1,:),SYMOUT_MRC(1,:));


% Using MRC and CMA or LMS
% mu      = [1e-3 1e-4 1e-6];
% ntaps   = [13 13 13];
% errid   = [1 1 1];
% iter    = [30 0 0];
% 
% polm    = 0;
% applms  = 1;
% CMAOUT  = DspAlg.PolarizationDemux(DSPIN,mn,sps,polm,mu,ntaps,errid,iter,applms);
% scatterplot(CMAOUT); title('w/ MRC');
% 
% polm    = 1;
% applms  = 1;
% CMAOUT  = DspAlg.PolarizationDemux(DSPIN,mn,sps,polm,mu,ntaps,errid,iter,applms);
% scatterplot(CMAOUT(:,1)); title('w/o MRC');
% % scatterplot(CMAOUT(:,2));

% [er,ec] = matlab_bert(CPEOUT,mn,OTX.BitReference,'gray')
% [er,ec] = matlab_bert(CPEOUT,mn,OTX.BitReference,'binary')
% [er,ec] = matlab_bert(CPEOUT,mn,OTX.BitReference,'diff')

