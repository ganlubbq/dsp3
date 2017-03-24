clear

osnr = 1 : 0.5 : 13;

for ii = 1 : length(osnr)
    vSet.nFrm = 1000;
    vSet.Hyd90 = 90;        % degree
    vSet.ADCfs = 2;         % samples per symbol
    vSet.DSPmode = 0;       % 0-offline; 1-real time
    vSet.DSPmemLen = 200;   % number of frames
    vSet.linewidth = 500e3;
    vSet.osnr = osnr(ii);
    vSet.baudrate = 30e9;
    vSet.bitpersym = 2;
    vSet.modFormat = 'QPSK';
    vSet.freqOffset = 1.0e9;
    vSet.psFiltType = 'Nyquist';
    %%%%%%%%%%%%%%%%%%%%%
    vM = go(vSet);
    %%%%%%%%%%%%%%%%%%%%%
    ber(ii) = vM.BER;
    snr(ii) = vM.SNR;
end

TheoreticalBERv2(2);

hold on;
semilogy(snr, ber, 'rs-.');