%HIST2 Plot color histogram and meshgrid of 2D complexed signal,
% return the handle of histogram figure and local peak of it, cmap defines
% the color map.
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

function [hs local_max] = Hist2(x1,cmap)
%
% check the input parameters
%
if nargin < 2
    cmap = 'hot';
end

%
% define grid
%
lim = 1.8;
coar_edge = linspace(-lim,lim,80);
fine_edge = linspace(-lim,lim,800);

%
% slight over normalization
%
x1 = DspAlg.Normalize(x1,16) / 3;

%
% generate histogram
%
N = length(coar_edge);
z1 = zeros(N,N);
[~,cbin] = histc(real(x1),coar_edge);
[~,dbin] = histc(imag(x1),coar_edge);
cbin = cbin + 1;
dbin = dbin + 1;
for ii = 1:length(x1)
    z1(cbin(ii),dbin(ii)) = z1(cbin(ii),dbin(ii)) + 1;
end

%
% find local peaks
%
[xc,yc] = meshgrid(coar_edge);
[xf,yf] = meshgrid(fine_edge);
z2 = interp2(xc,yc,z1,xf,yf);
z3 = smooth2a(z2,30,30);
xy = FastPeakFind(z3);
xy = (xy-400)/400*lim - 0.015;
local_max = xy(:,1) + 1i*xy(:,2);

%
% plot local peaks
%
hs = scatterplot(local_max);

%
% plot 2D color histogram
%
surface(coar_edge,coar_edge,z1);
axis([-lim lim -lim lim]);
shading interp

%
% setup color map
%
map = colormap(cmap);
% map(1,:) = [0 0 0];
set(hs,'Colormap',map);
set(get(hs,'CurrentAxes'),'Color',map(1,:))

%
% plot 3D mesh grid
%
hm = figure;
mesh(xf,yf,z3);
view(-20,60);
set(hm,'Colormap',map);

%
% manage figure positions
%
h2_p = get(hs,'Position');
h3_p = get(hm,'Position');
h2_p(1) = h2_p(1) - h2_p(3)/2 - 5;
h3_p(1) = h3_p(1) + h2_p(3)/2 + 5;
set(hs,'Position',h2_p);
set(hm,'Position',h3_p);

return

