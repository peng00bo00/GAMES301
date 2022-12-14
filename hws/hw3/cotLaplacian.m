function L = cotLaplacian(V, F)
%% Cotangent-Laplacian operator

nV = size(V, 1);

E = reshape(V(F(:, [2, 3, 1]), :) - V(F, :), [size(F), 3]);
cots = -dot(E(:, [2, 3, 1], :), E(:, [3, 1, 2], :), 3) ./ vecnorm(cross(E(:, [2, 3, 1], :), E(:, [3, 1, 2], :), 3), 2, 3);

Gcot = sparse(F, F(:, [2, 3, 1]), cots, nV, nV);
L = Gcot + Gcot.';

L = diag(sparse(sum(L, 1))) - L;

end