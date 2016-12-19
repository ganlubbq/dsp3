function bit = slicerBPSK(sym)
% Convert square mQAM symbols to bits in rows according to gray mapping.
% The decimal order of uncoded symbol in contellation is from topleft to
% bottomright by columns
% 
% Example: bit = slicerBPSK(sym)
% 
% Input: 
%       sym         - input symbols in column
% 
% Reference: 
% 
% Note: 
% 
% See Also: 

% put dec through mapper than do dec2bin
bit = real(sym)>0;

return

