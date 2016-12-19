function [bb bg] = MakeDecision(x,mn)
%MAKEDECISION Ouput binary or gray-decoded bits for qpsk and 16qam
%
%   See Also: sym2bit,bit2sym,CoderDriver
%
%   symbol map for 16-qam (binary encode):
%       [4 8 12 16; ...
%        3 7 11 15; ...
%        2 6 10 14; ...
%        1 5 9 13];
%
%   symbol map for 16-qam (gray encode):
%       [3 7 11 15; ...
%        4 8 12 16; ...
%        2 6 10 14; ...
%        1 5 9 13];
%
%   symbol map for qpsk:
%       [2 4;...
%        1 3]

%   Copyright 2012

switch mn
    case 2
        rx = sign(real(x));
    case 4
        rx = [sign(real(x)), sign(imag(x))];
    case 16
        bound = 2;
        rx = [sign(real(x)), sign(real(x)-sign(real(x))*bound),...
              sign(imag(x)), sign(imag(x)-sign(imag(x))*bound)];
    case 64
        bound = [2,4,6];
        t1 = [sign(real(x)), ...
            sign(real(x)-sign(real(x))*bound(1)), ...
            sign(real(x)-sign(real(x))*bound(2)), ...
            sign(real(x)-sign(real(x))*bound(3))];
        t2 = [sign(imag(x)), ...
            sign(imag(x)-sign(imag(x))*bound(1)), ...
            sign(imag(x)-sign(imag(x))*bound(2)), ...
            sign(imag(x)-sign(imag(x))*bound(3))];
        rx = [t1 t2];
end
bb = sign(rx+1);
bg = DspAlg.sym2bit(DspAlg.bit2sym(bb,mn),mn);

