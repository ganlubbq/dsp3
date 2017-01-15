%% TEST SCRIPT FOR PARAMETERIZED PDF OF CLASSICAL ESTIMATION
% SLICE ALONG FIXED DATA WOULD GENERATE A LIKELIHOOD OF UNKNOWN PARAMTER
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
%% Unknown parameter is deterministic
clear

% unknown parameter
theta = -1:0.05:1;

% likelihood function
PDF = 1/sqrt(2*pi*0.1)*exp(-theta.^2/(2*0.1));
PDF = PDF';

% varying with data
PDF3 = [];
for jj = 1:length(theta)
    PDF3 = [PDF3 PDF+0.2*rand(length(theta),1)];
end

% nominal data
x = -1:0.05:1;

figure;
mesh(theta,x,PDF3); view([-134 28]);
xlabel('x'); ylabel('\theta'); zlabel('PDF');
