% DESCRIPTION
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

function mf_name = getModemFormat( mn )

switch mn
    case 1
        mf_name = 'ASK';
    case 2
        mf_name = 'BPSK';
    case 4
        mf_name = 'QPSK';
    case 8
        mf_name = '8PSK';
    case 16
        mf_name = '16-QAM';
    case 32
        mf_name = '32-QAM';
    case 64
        mf_name = '64-QAM';
    case 128
        mf_name = '128-QAM';
    case 256
        mf_name = '256-QAM';
    case 512
        mf_name = '512-QAM';
    case 1024
        mf_name = '1024-QAM';
end

return

