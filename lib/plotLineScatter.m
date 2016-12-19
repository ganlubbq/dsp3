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

function h = plotLineScatter(x,y)

group = {'bs-','gd-','ro-','c>-','m<-','kx-',...
    'bs--','gd--','ro--','c>--','m<--','kx--'};

h = figure; hold on; grid on;

if nargin == 2
    for ii = 1:size(y,2)
        plot(x(:),y(:,ii),group{ii},'Linewidth',2,'MarkerSize',12);
    end
    legend toggle
end

if nargin == 1
    for ii = 1:size(x,2)
        plot(x(:,ii),group{ii},'Linewidth',2,'MarkerSize',12);
    end
    legend toggle
end

hold off

return

