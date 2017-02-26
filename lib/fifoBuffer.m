function [newstate] = fifoBuffer(state,data)
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
% Copyright default

% number of columns, data length
N = size(data,2);

if N > size(state,2)
    warning('data length larger than buffer length, lossing data');
    N = size(state,2);
end

% fifo
newstate = [state(:,N+1:end) data(:,1:N)];

return
