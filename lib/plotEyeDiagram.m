% DESCRIPTION
% 
% Example: h = plotEyeDiagram(x, baudRate, samplingRate, xMode)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function h = plotEyeDiagram(h, x, baudRate, samplingRate, xMode)

if length(x) <= 1000
    N = length(x);
else
    N = 1000;
end
numEyes = 2;
period = numEyes;
offset = 0;
plotstr = 'b';
sps = samplingRate/baudRate;
% need 120 sps for smoothing eye diagram
upfactor = round(120/sps);
spsFinal = upfactor*sps;

% h = commscope.eyediagram(...
%     'SamplingFrequency', samplingRate, ...
%     'SamplesPerSymbol', sps * upfactor, ...
%     'MinimumAmplitude', 0, ...
%     'MaximumAmplitude', 1, ...
%     'PlotType', '2D Color', ...
%     'ColorScale', 'log', ...
%     'RefreshPlot', 'on'...
%     );

switch xMode
    case 'optical'
        if isreal(x)
            xd = resample(x(1:N),upfactor,1,2); % 2nd order filter, smoother than default 10th order
            data = (abs(xd)./max(abs(xd))).^2;
        else
            xi = real(x);
            xq = imag(x);
            xid = resample(xi(1:N),upfactor,1,2);
            xqd = resample(xq(1:N),upfactor,1,2);
            data = xid + 1j*xqd;
            data = (abs(data)./max(abs(data))).^2;
        end
        h = eyediagram(data,spsFinal*numEyes,period,offset,plotstr,h);
    case 'electrical'
        if isreal(x)
            xd = resample(x(1:N),upfactor,1,2);
            data = xd./max(abs(xd));
        else
            xi = real(x);
            xq = imag(x);
            xid = resample(xi(1:N),upfactor,1,2);
            xqd = resample(xq(1:N),upfactor,1,2);
            data = xid + 1j*xqd;
            data = abs(data)./max(abs(data));
        end
        h = eyediagram(data,spsFinal*numEyes,period,offset,plotstr,h);
    otherwise
        error('input mode not supported');
end

return

