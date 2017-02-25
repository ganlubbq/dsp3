function [U] = orthGramSchmidt(V)
%ORTHGRAMSCHMITDT implements the stabilized Gram–Schmidt orthonormalization. The
%columns of matrix V are replaced by orthonormal vectors (columns of U)
%which span the same subspace.
%
% Classic Gram-Schmidt process:
% uk = vk - proj_of_vk_on_u1 - proj_of_vk_on_u2 - ... -
% proj_of_vk_on_u(k-1)
%
% Modified (stablized) Gram-Schmidt process computes each perpendicular
% component iteratively so that each perpendicular component is orthogonal
% to any error introduced in each step:
% 1st: uk = vk - proj_of_vk_on_u1;
% 2nd: uk = uk - proj_of_uk_on_u2;
% ....

n = size(V,1);
k = size(V,2);
U = zeros(n,k);
U(:,1) = V(:,1) / sqrt(V(:,1)' * V(:,1));
for i = 2 : k
    U(:,i) = V(:,i);
    for j = 1:i-1
        U(:,i) = U(:,i) - (U(:,i)' * U(:,j)) * U(:,j);
    end
    U(:,i) = U(:,i) / sqrt(U(:,i)' * U(:,i));
end

end