function L = cotLaplacian(V, F)
%% Build (sparse) cotangent-Laplacian matrix
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      L[nV, nV]: sparse cotangent-Laplacian matrix

nV = size(V, 1);

%% edges
Es = reshape(V(F(:, [2, 3, 1]), :) - V(F, :), [size(F), 3]);

%% trigonometry
coss =-dot(Es(:, [2, 3, 1], :), Es(:, [3, 1, 2], :), 3);
sins = vecnorm(cross(Es(:, [2, 3, 1], :), Es(:, [3, 1, 2], :), 3), 2, 3);
cots = coss ./ sins;

%% cotangent-Laplacian adjacency
L = sparse(F, F(:, [2, 3, 1]), cots, nV, nV);
L = 0.5*(L + L');

%% L = D - A
L = diag(sparse(sum(L, 1))) - L;

end