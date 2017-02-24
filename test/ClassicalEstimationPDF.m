%% TEST SCRIPT FOR PARAMETERIZED PDF OF CLASSICAL ESTIMATION
% SLICE ALONG FIXED DATA WOULD GENERATE A LIKELIHOOD FUNCTION OF UNKNOWN PARAMTER
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
%% Unknown parameter is deterministic
clear

% nominal data
x = -1:0.05:1;

% the deterministic unknown parameter, gaussian mean
theta = -1:0.05:1;

% % 
% PDF = 1/sqrt(2*pi*0.1)*exp(-(x.-theta).^2/(2*0.1));
% PDF = PDF';

% varying with data
PDF3 = [];
for jj = 1:length(theta)
    PDF = 1/sqrt(2*pi*0.1)*exp(-(x-theta(jj)).^2/(2*0.1));
    PDF = PDF';
    PDF3 = [PDF3 PDF];
end

% plot the parameterized pdf
h1=figure;
mesh(theta,x,PDF3); view([-64 19]);
xlabel('x'); ylabel('\theta'); zlabel('PDF');

% the likelihood function can be obtained by fix data x and varying unknown
% parameter
h2=figure;
plot(theta,PDF3(10,:)); grid on

mngFigureWindow(h1,h2);

