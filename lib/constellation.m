function c = constellation(mn)
%CONSTELLATION Give the canonical form of uncoded constellations for mn-QAM
%

h = modem.qammod('M',mn);

cr = h.Constellation;

if mn == 8
    c = [-2, -2j, 2j, 2, 3+3j, 3-3j, -3-3j, -3+3j];
else
    c = cr(:);
end

return


