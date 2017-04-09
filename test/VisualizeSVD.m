% Visualize sigular value decomposition of matrix. For an arbitrary matrix,
% one can always find a set of orthonormal vectors V which will be
% transformed to another set of orthonormal vectors U. The scaling factor
% for each vector is called the singular values.
clear

% a sheer matrix
A = [1, 1; 0, 1];

% input space v is a circle
x = -1:0.01:1;
y = sqrt(1 - x.^2);
v = [x(:), y(:); flipud(x(:)), flipud(-y(:))];

% linear transform
u = A * v.';
u = u.';

% visualize
linewidth = 2;
h = figure('Position', [483 188 572 539]);
hold on
grid on
box on
axis([-1.5 1.5 -1.5 1.5]);
plot(v(:,1), v(:,2), 'b-', 'LineWidth', linewidth);
plot(u(:,1), u(:,2), 'g-.', 'LineWidth', linewidth);

% svd
[U, S, V] = svd(A);

% visualize the orthonormal basis
line([0, V(1,1)], [0, V(2,1)], 'LineWidth', linewidth);
line([0, V(1,2)], [0, V(2,2)], 'LineWidth', linewidth);
line([0, U(1,1) * S(1,1)], [0, U(2,1) * S(1,1)], 'color', 'g', 'LineWidth', linewidth);
line([0, U(1,2) * S(2,2)], [0, U(2,2) * S(2,2)], 'color', 'g', 'LineWidth', linewidth);
