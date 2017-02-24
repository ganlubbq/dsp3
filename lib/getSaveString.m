function [varargout] = getSaveString(Des, patExtdir)
%GETSAVESTRING Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    patExtdir = '';
end
if nargin < 1
    Des = '';
end
if isempty(Des)
    Des = 'missing description';
end

[DateMarker, TimeMarker, Logtime] = getClock();

if ~isempty(patExtdir)
    if ~strcmp(patExtdir(end),'\')
        patExtdir = [patExtdir,'\'];
    end
end

Savedir = [patExtdir,DateMarker,'\'];

if ~exist(Savedir,'dir')
    mkdir(Savedir);
end

Logfile = 'log.txt';
Datfile = [TimeMarker,'.mat'];
Descrip = Des;

varargout{1} = Savedir;
varargout{2} = Logfile;
varargout{3} = Datfile;
varargout{4} = Descrip;
varargout{5} = Logtime;

return
