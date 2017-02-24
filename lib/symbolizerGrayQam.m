function sym = symbolizerGrayQam(bit)
% Convert bits in rows to square mQAM symbols according to gray mapping.
% The decimal order of uncoded symbol in contellation is from topleft to
% bottomright by columns
% 
% Example: sym = symbolizerGrayQam(bit)
% 
% Input: bit - input bits in rows, the number of rows reps number of bit
% per symbol
% 
% Reference: 
% 
% Note: Gray mapping only
% 
% See Also: 

% size of alphabet
M = 2 ^ size(bit,1);

nSample = size(bit,2);

% topleft to bottomright by columns
c = constellation(M);

% special case for bpsk
if M == 2
    sym = c(bit+1); return
end

% convert bit in rows to dec by LSB to MSB order
twos = 2.^(0:log2(M)-1);
twos = repmat(twos,nSample,1).';
dec = sum(bit.*twos)+1;

% integer mapper of gray coded m-qam
mapper = mapint(M);

% put dec through int mapper then do constellation mapping
sym = c(mapper(dec));

return

% gray integer mapper for m-qam
function order = mapint(mn)

switch mn
    case 4
        order = [1,2,3,4;];
    case 16
        order = [1,2,4,3,5,6,8,7,13,14,16,15,9,10,12,11;];
    case 64
        order = [1,2,4,3,8,7,5,6,9,10,12,11,16,15,13,14,25,26,28,27,32,31,29,30,17,18,20,19,24,23,21,22,57,58,60,59,64,63,61,62,49,50,52,51,56,55,53,54,33,34,36,35,40,39,37,38,41,42,44,43,48,47,45,46;];
    case 256
        order = [1,2,4,3,8,7,5,6,16,15,13,14,9,10,12,11,17,18,20,19,24,23,21,22,32,31,29,30,25,26,28,27,49,50,52,51,56,55,53,54,64,63,61,62,57,58,60,59,33,34,36,35,40,39,37,38,48,47,45,46,41,42,44,43,113,114,116,115,120,119,117,118,128,127,125,126,121,122,124,123,97,98,100,99,104,103,101,102,112,111,109,110,105,106,108,107,65,66,68,67,72,71,69,70,80,79,77,78,73,74,76,75,81,82,84,83,88,87,85,86,96,95,93,94,89,90,92,91,241,242,244,243,248,247,245,246,256,255,253,254,249,250,252,251,225,226,228,227,232,231,229,230,240,239,237,238,233,234,236,235,193,194,196,195,200,199,197,198,208,207,205,206,201,202,204,203,209,210,212,211,216,215,213,214,224,223,221,222,217,218,220,219,129,130,132,131,136,135,133,134,144,143,141,142,137,138,140,139,145,146,148,147,152,151,149,150,160,159,157,158,153,154,156,155,177,178,180,179,184,183,181,182,192,191,189,190,185,186,188,187,161,162,164,163,168,167,165,166,176,175,173,174,169,170,172,171;];
    otherwise
        warning('unsupported format');
        keyboard;
end

return