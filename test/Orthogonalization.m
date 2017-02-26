%% SCRIPT FOR TESTING VECTOR ORTHOGONALIZATION USING GRAM-SCHMIDT PROCESS
%

clear

xi = rand(10000,1)-0.5;
xq = rand(10000,1)-0.5;
% figure; plot(xi+1i*xq,'.'); grid on

x = [xi xq];
x = x.';

% rotate x by 15 degree
V = [1 cos(75/180*pi); 0 sin(75/180*pi)];
r = V * x;

% verify the orthogonality
inner_product = r(1,:) * r(2,:)'

h1 = figure; plot(r(1,:),r(2,:),'.'); grid on
xlim([-1 1]); ylim([-1 1]);

%
u = orthGramSchmidt(r.');
h2 = figure; plot(u(:,1),u(:,2),'.'); grid on

% verify the orthogonality
inner_product = u(:,1)' * u(:,2)

mngFigureWindow(h1, h2);

% EOF