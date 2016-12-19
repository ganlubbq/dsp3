% FIFO buffering data in rows
%
% Example: [newstate] = fifoBuffer(state,data)
% 
% Input: 
%       state       - current buffer
%       data        - for buffering
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

function [newstate] = fifoBuffer(state,data)

N = size(data,2);

if N > size(state,2)
    warning('data length larger than buffer length');
    N = size(state,2);
end

data1 = state(:,N+1:end);

newstate = [data1 data(:,1:N)];

return


