% Convert bits in rows to BPSK symbols according to gray mapping.
% The decimal order of uncoded symbol in contellation is from topleft to
% bottomright by columns
% 
% Example: sym = symbolizerBPSK(bit)
% 
% Input: 
%       bit - input bits in rows, the number of rows reps number of bit per
%             symbol
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright default

function sym = symbolizerBPSK(bit)

% topleft to bottomright by columns
c = constellation(2);

sym = c(bit+1);

return

% canonical uncoded constellation
function c = constellation(mn)

h = modem.qammod('M',mn);

cr = h.Constellation;

c = cr(:);

return

return

