% Convert bits in rows to symbols according to certain mapping rules. The
% mapping order is from topleft to bottomright by columns
% 
% Example: 
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

function sym = symbolizer(bit,mapint)

% size of alphabet
M = size(bit,1);

nSample = size(bit,2);

% topleft to bottomright by columns
c = constellation(M);

% convert bit in rows to int
twos = 2.^(0:M-1);
twos = repmat(twos,nSample,[]).';
int = sum(bit.*twos);

% mapping
sym = c(int);

return

