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
% Copyright 2015 Default

function h = plotEyeDiagram(x, period, xMode)

if ~iscolumn(x)
    warning('input has to be a column vector'); keyboard;
end

offset = 0;

% h = commscope.eyediagram(...
%     'SamplingFrequency', samplingRate, ...
%     'SamplesPerSymbol', sps * upfactor, ...
%     'MinimumAmplitude', 0, ...
%     'MaximumAmplitude', 1, ...
%     'PlotType', '2D Color', ...
%     'ColorScale', 'log', ...
%     'RefreshPlot', 'on'...
%     );

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

