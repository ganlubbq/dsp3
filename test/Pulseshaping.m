%% TEST SCRIPT FOR TESTING VARIOUS TYPES OF PULSE-SHAPING FILTERING
% Compare raised cosine pulse-shaping filter in both frequency and time
% domain
%%

nSymbol = 2^16;

% NUMBER OF BIT PER SYMBOL
k = 2;

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
signal_power = sum(abs(sym).^2)/nSymbol

sps = 16;

nSamples = sps*nSymbol;

% Upsampling
if ~iscolumn(sym)
    sym = sym.';
end
temp = sym * ones(1,sps);
temp(:,2:end) = 0; % insert zeros
temp = temp.';
% temp = repmat(sym,sps,1);
sym_upsampled = temp(:);

% get a freq domain raised cosine filter response
Rs = 1;
Fs = sps;
freqVect = getFFTGrid(nSamples,Fs);
alpha = 0.35; mode = 0;
H = calcRcosResponse(nSamples,Fs,Rs,alpha,mode);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)).*H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)).*H));

sym_filtered = sym_upsampled_i + 1j*sym_upsampled_q;

% display signal power after pulse shaping
signal_power_rcos_freq = sum(abs(sym_filtered).^2)/nSamples

h1 = plotEyeDiagram(sym_filtered(1:8192),2*sps,'e');

% design a rcos digital filter
span = 10;
shape = 'normal';
h = rcosdesign(alpha,span,sps,shape);
delay = span*sps/2;

if ~iscolumn(sym)
    sym = sym.';
end
sym_upsampled_i = upfirdn(real(sym), h, sps);
sym_upsampled_q = upfirdn(imag(sym), h, sps);

% if iscolumn(sym_upsampled)
%     sym_upsampled = sym_upsampled.';
% end
% sym_upsampled_i = firfilt(real(sym_upsampled), h, 'overlap-save','same');
% sym_upsampled_q = firfilt(imag(sym_upsampled), h, 'overlap-save','same');

sym_filtered = sym_upsampled_i + 1j*sym_upsampled_q;

% display signal power after pulse shaping
signal_power_rcos_time = sum(abs(sym_filtered).^2)/nSamples

if ~iscolumn(sym_filtered)
    sym_filtered = sym_filtered.';
end
h2 = plotEyeDiagram(sym_filtered((1:8192)+delay),2*sps,'e');

mngFigureWindow(h1,h2);

