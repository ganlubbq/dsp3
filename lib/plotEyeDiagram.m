function h = plotEyeDiagram(x, period, xMode)
% DESCRIPTION
% 
% Example: h = plotEyeDiagram(x, baudRate, samplingRate, xMode)
% 
% Input: period - number of samples per eye
% 

if ~iscolumn(x)
    error('input has to be a column vector');
end

% fancy one...
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
xds = x(1 : end - tmp);
xrs = reshape(xds, period, []);

switch xMode
    case 'o'
        h = figure; plot(0 : period - 1, abs(xrs).^2, 'b'); grid on
        xlim([0, period]);
    case 'e'
        if isreal(x)
            h = figure; plot(0 : period - 1, xrs, 'b'); grid on
            xlim([0, period]);
        else
            h = figure;
            subplot(211); plot(0 : period - 1, real(xrs), 'b'); grid on; xlim([0, period]);
            subplot(212); plot(0 : period - 1, imag(xrs), 'g'); grid on; xlim([0, period]);
        end
    otherwise
        warning('input data mode not supported'); keyboard;
end

return
