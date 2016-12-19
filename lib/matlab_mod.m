function [bint,h] = matlab_mod(x,mn)
%
% # input/bit | binary output
% # input/integer | symbol output
%
% # symorder/binary
% # symorder/gray
% # symorder/user-defined
%

h = modem.qammod;
h.M = mn;
h.PhaseOffset = 0;
h.SymbolOrder = 'binary';
h.InputType = 'integer';

bint = h.modulate(x);
