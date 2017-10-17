function bit = slicer_mqam(sym, M)
% Convert square mQAM symbols to bits in rows according to gray mapping.
% The decimal order of uncoded symbol in contellation is from topleft to
% bottomright by columns
% 
% Input: 
%       bit         - output bits in rows
%       M           - alphabet size
%       sym         - input symbols in column
% 
% Reference: 
% 
% h = modem.qamdemod;
% h.M = mn;
% h.PhaseOffset = 0;
% h.SymbolOrder = 'binary';
% 
% h.OutputType = 'integer';
% bint = h.demodulate(x);
% b = h.Constellation(bint+1);
% 
% Note: 
% 
% See Also: symbolizerGrayQam

% special case of bpsk
if M == 2
    bit = real(sym) > 0; return
end

% number of rows of output bit
k = log2(M);

% normalize symbols
sym = normalizeQam(sym, M);

% convert symbol to decimal index from topleft to bottomright by columns,
% i.e. sym->[I,~Q]
decndx = slicer(sym, M);

% integer mapper of gray coded m-qam
mapper = mapint(M);

% put dec through mapper than do dec2bin
bit = dec2bin(mapper(decndx), k);

return

function order = mapint(mn)

switch mn
    case 4
        order = [0,1,2,3;];
    case 16
        order = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
    case 64
        order = [0,1,3,2,6,7,5,4,8,9,11,10,14,15,13,12,24,25,27,26,30,...
            31,29,28,16,17,19,18,22,23,21,20,48,49,51,50,54,55,53,52,...
            56,57,59,58,62,63,61,60,40,41,43,42,46,47,45,44,32,33,35,...
            34,38,39,37,36;];
    case 256
        order = [0,1,3,2,6,7,5,4,12,13,15,14,10,11,9,8,16,17,19,18,22,...
            23,21,20,28,29,31,30,26,27,25,24,48,49,51,50,54,55,53,52,...
            60,61,63,62,58,59,57,56,32,33,35,34,38,39,37,36,44,45,47,...
            46,42,43,41,40,96,97,99,98,102,103,101,100,108,109,111,110,...
            106,107,105,104,112,113,115,114,118,119,117,116,124,125,127,...
            126,122,123,121,120,80,81,83,82,86,87,85,84,92,93,95,94,90,...
            91,89,88,64,65,67,66,70,71,69,68,76,77,79,78,74,75,73,72,192,...
            193,195,194,198,199,197,196,204,205,207,206,202,203,201,200,...
            208,209,211,210,214,215,213,212,220,221,223,222,218,219,217,...
            216,240,241,243,242,246,247,245,244,252,253,255,254,250,251,...
            249,248,224,225,227,226,230,231,229,228,236,237,239,238,234,...
            235,233,232,160,161,163,162,166,167,165,164,172,173,175,174,...
            170,171,169,168,176,177,179,178,182,183,181,180,188,189,191,...
            190,186,187,185,184,144,145,147,146,150,151,149,148,156,157,...
            159,158,154,155,153,152,128,129,131,130,134,135,133,132,140,...
            141,143,142,138,139,137,136;];
    otherwise
        error('unsupported format');
end

return

function [ndx] = slicer(x,mn)

x = x(:);

switch mn % MSB first
    
    case 4
        bit = [msign(real(x)), ~msign(imag(x))];
        
    case 16
        bound = 2;
        bit = [msign(real(x)), msign(real(x)-ksign(real(x))*bound),...
              ~msign(imag(x)), ~msign(imag(x)-ksign(imag(x))*bound)];
          
    case 64 % 64qam slicing is combination of qpsk and 16qam, i.e.
    % [4qam 16qam 16qam 4qam 16qam 16qam]
        bound = 2;
        
        % 4qam header
        bit1 = [msign(real(x)), ~msign(imag(x))];
        
        % convert 64qam symbol in one quater to 16qam
        x = x-(ksign(real(x))*4+ksign(imag(x))*4j);
        
        % slicing 16qam
        bit2 = [msign(real(x)), msign(real(x)-ksign(real(x))*bound),...
            ~msign(imag(x)), ~msign(imag(x)-ksign(imag(x))*bound)];
        
        bit = [bit1(:,1) bit2(:,1:2) bit1(:,2) bit2(:,3:4)];
    
    case 256
        % TBA
    otherwise
        error('unsupported format');
end

% convert bit in rows to dec by MSB to LSB order
twos = 2.^(log2(mn)-1 : -1 : 0);
twos = repmat(twos, length(x), 1);
ndx = sum(bit .* twos, 2) + 1;

return

% return 1 if x>0 and 0 if x<=0
function s = msign(x)

s = (x > 0);

return

% return 1 if x>0 and -1 if x<=0
function s = ksign(x)

s = (x > 0) * 2 - 1;

return

% convert decimal index to bits in k rows by LSB to MSB order
function bit = dec2bin(dec,k)

dec = dec(:).';

bit = zeros(k, length(dec));

for ii=1:k
    bit(ii,:) = mod(dec, 2);
    dec = floor(dec / 2);
end

return
