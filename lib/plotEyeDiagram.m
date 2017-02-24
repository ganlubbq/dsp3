function h = plotEyeDiagram(x, period, xMode)
% DESCRIPTION
% 
% Example: h = plotEyeDiagram(x, baudRate, samplingRate, xMode)
% 
% Input: period - number of samples per eye
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 

if ~iscolumn(x)
    error('input has to be a column vector');
end

offset = 0;

% fancy one
% h = commscope.eyediagram(...
%     'SamplingFrequency', samplingRate, ...
%     'SamplesPerSymbol', sps * upfactor, ...
%     'MinimumAmplitude', 0, ...
%     'MaximumAmplitude', 1, ...
%     'PlotType', '2D Color', ...
%     'ColorScale', 'log', ...
%     'RefreshPlot', 'on'...
%     );

% get rid of extra data
tmp = mod(length(x), period);
x = x(1:end-tmp);

x = reshape(x,period,[]);

switch xMode
    case 'o'
        h = figure;
        plot(abs(x).^2,'b'); grid on
    case 'e'
        if isreal(x)
            h = figure;
            plot(x,'b'); grid on
        else
            h = figure;
            subplot(211);plot(real(x),'b');grid on
            subplot(212);plot(imag(x),'g');grid on
        end
    otherwise
        warning('input data mode not supported'); keyboard;
end

return
