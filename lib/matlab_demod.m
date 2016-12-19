function [bbit,bint,h] = matlab_demod(x,mn,method)
%
%   input/symbol
%   
%   slicer + sym2bit/sym2de
%

if isreal(x)
    error('real data is not supported.')
end

x = DspAlg.Normalize(x,mn);

h = modem.qamdemod;
h.M = mn;
h.PhaseOffset = 0;
h.SymbolOrder = method;

h.OutputType = 'integer';
bint = h.demodulate(x);

h.OutputType = 'bit';
bbit = h.demodulate(x);