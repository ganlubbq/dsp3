function [h] = polarizationAnalyzer(ex,ey,varargin)
% Polarization analyzer plot SOP as a function of frequency on a Poincare
% sphere. At any single frequency, the signal is completely polarized.
% 
% Example: [h] = polarization_analyzer(ex,ey,varargin)
% 
% Input: 
%       ex - must be in row
%       ey - must be in row
% 
% Reference: 
% 
% Note: Input field can be either time or freq domain
% 
% See Also: spectrumAnalyzer

s = j2s(ex,ey);

h = poincare_sphere();
plot3(s(1,:),s(2,:),s(3,:),varargin{:});

return

function s_out = j2s(ex,ey)
% Transfer Jones vector to normalized Stokes vector

% input Jones vector either in time or freq domain
e_in = [ex(:) ey(:)].';

% Pauli matrix
sigma = {[1,0;0,-1],[0,1;1,0],[0,-1i;1i,0]};

data_length = length(e_in(1,:));
s_out = zeros(3,data_length);

for loop_1 = 1:3
    e_tmp = sigma{loop_1} * e_in;
    for loop_2 = 1:data_length
        s_out(loop_1,loop_2) = e_in(:,loop_2)' * e_tmp(:,loop_2);
    end
end

% normalize
s0 = sqrt(s_out(1,:).^2 + s_out(2,:).^2 + s_out(3,:).^2);
s_out = [s_out(1,:)./s0;s_out(2,:)./s0;s_out(3,:)./s0];

return

function h = poincare_sphere()
% Plot a Poincar sphere

leftbottom = 80;
width = 800;
height = 600;

h = figure('visible','on');

set(gcf,'position', ...
    [leftbottom,leftbottom,leftbottom+width,leftbottom+height])

sphere
axis('equal','off')
colormap(gray)
shading interp
alpha(0.8)
view([135,25])

hold on

plot3([0,2],[0,0],[0,0],'w-.','linewidth',.5)
plot3([0,0],[0,2],[0,0],'w-.','linewidth',.5)
plot3([0,0],[0,0],[0,1.5],'w-.','linewidth',.5)

theta_tmp = linspace(0,2*pi,60);
x_tmp = cos(theta_tmp);
y_tmp = sin(theta_tmp);
z_tmp = zeros(1,length(theta_tmp));
plot3(x_tmp,y_tmp,z_tmp,'w-','linewidth',.5)

x_tmp = cos(theta_tmp);
y_tmp = zeros(1,length(theta_tmp));
z_tmp = sin(theta_tmp);
plot3(x_tmp,y_tmp,z_tmp,'w-','linewidth',.5)

x_tmp = zeros(1,length(theta_tmp));
y_tmp = cos(theta_tmp);
z_tmp = sin(theta_tmp);
plot3(x_tmp,y_tmp,z_tmp,'w-','linewidth',.5)

clear *tmp

text(2.1,0,0,'S_1','FontWeight','b','FontSize',14)
text(0,2.1,0,'S_2','FontWeight','b','FontSize',14)
text(0,0,1.6,'S_3','FontWeight','b','FontSize',14)

return
