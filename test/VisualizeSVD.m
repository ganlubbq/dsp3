
clear

% a sheer matrix
A = [1 1; 0 1];

% input space v is a circle
x = -1:0.01:1;
y = sqrt(1 - x.^2);
v = [x(:), y(:); flipud(x(:)), flipud(-y(:))];

% linear transform
u = A * v.';
u = u.';

% visualize
linewidth = 2;
figure;
hold on
grid on
box on
axis([-1.5 1.5 -1.5 1.5]);
axis square
plot(v(:,1), v(:,2), 'b-', 'LineWidth', linewidth);
plot(u(:,1), u(:,2), 'g-.', 'LineWidth', linewidth);

% svd
[U, S, V] = svd(A);

% visualize the orthonormal basis
line([0, V(1,1)], [0, V(2,1)], 'LineWidth', linewidth);
line([0, V(1,2)], [0, V(2,2)], 'LineWidth', linewidth);
line([0, U(1,1) * S(1,1)], [0, U(2,1) * S(1,1)], 'color', 'g', 'LineWidth', linewidth);
line([0, U(1,2) * S(2,2)], [0, U(2,2) * S(2,2)], 'color', 'g', 'LineWidth', linewidth);
