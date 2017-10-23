function [newstate] = fifo_buffer(state, data)
% FIFO buffering data in rows
% 
% [newstate] = fifo_buffer(state, data)
%	state : current buffer
%   data  : for buffering

% number of columns, data length
ncol = size(data, 2);
if ncol > size(state,2)
    warning('data length larger than buffer length, lossing data');
    ncol = size(state,2);
end

% fifo
newstate = [state(:, ncol+1:end), data(:, 1:ncol)];
return
