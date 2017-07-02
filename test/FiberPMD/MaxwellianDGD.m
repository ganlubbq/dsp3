% This script draws the maxwellian distribution of the fiber DGD, which is
% defined as the length of the PMD vector.
%
% Magnus Karlsson, "Probability density functions of the DGD in optical
% fiber communication systems" JLT 2001
%
clear

% fiber length
L = 100e3;

% number of sections
N = 100;

% DGD parameter
DGD_Param = 5 * 1e-12 / sqrt(1000);

% mean DGD
DGD_Mean = DGD_Param * sqrt(L);

% DGD in each section
DGD_Section = DGD_Mean * sqrt(3 * pi /8 / N);

% Maxiwellian PDF
tau = 0 : 1e-12 : 150e-12;
pMxw = (3 * tau.^2 / (DGD_Section.^3 * N)) .* sqrt(6 / pi / N) .* exp(-3 * tau.^2 / (2 * N * DGD_Section.^2)) / 1e12;

% CDF
cMxw = cumsum(pMxw);

figure;
subplot(211); plot(tau, pMxw, 'b', 'LineWidth', 4); grid on
xlabel('DGD'); ylabel('Prob.');
legend('Maxwellian PDF');

subplot(212); plot(tau, cMxw, 'r', 'LineWidth', 4); grid on
xlabel('DGD'); ylabel('Prob.');
legend('Maxwellian CDF');
