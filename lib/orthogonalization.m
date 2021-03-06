function [y, theta] = orthogonalization(x)
%ORTHOGONALIZATION Gram-Schmidt orthogonalization procedure
% we have vector {r1,r2} and the angle between them is pi/2-theta, now,
% build two orthogonal vectors from {r1,r2}:
% =============================
% r1' = r1
% r2' = r2 - <r1r2>/<r1^2>*r1'
% =============================
for k = 1 : size(x, 2)
    r1 = real(x(:, k));
    r2 = imag(x(:, k));
    sin2theta = mean(r1 .* r2) / mean(r1.^2);
    r4 = (r2 - r1 * sin2theta) / sqrt(1 - sin2theta^2);
    y(:, k) = r1 + 1i * r4;
end
theta = asin(sin2theta) / pi * 180;
return
