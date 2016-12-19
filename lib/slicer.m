function b = slicer(x,mn)
%
%

h = modem.qamdemod;
h.M = mn;
h.PhaseOffset = 0;
h.SymbolOrder = 'binary';

h.OutputType = 'integer';
bint = h.demodulate(x);
b = h.Constellation(bint+1);