function [uv] = tutte(X, F, B)

% input: X (geometry), F (topology), B (boundary indices)
% output: parameterization for each vertex

nv = size(X, 1);       % number of vertices
nb = length(B);         % length of boundary
Gvv = sparse(F, F(:, [2, 3, 1]), 1, nv, nv); 
Gvv( Gvv | Gvv.' ) = 1;  % VV adjacent sparse matrix
D = diag(sum(Gvv));      % degree diagnal matrix (sparse)

Laplacian = D - Gvv;     % Laplacian of mesh

Xb = setBoundary(nb, 'circle');     % boundary condition
y = Gvv(:, B) * Xb;                 % right-hand side
inB = logical(zeros(1, nv));        % whether i is in boundary
inB(B) = 1;

M = Laplacian(~inB, ~inB);
b = y(~inB, :);                 % linear system exclude boundary

L = chol(M, 'lower');           % cholesky decomposition

uv = zeros(nv, 2);
uv(~inB, :) = [L' \ (L \ b(:, 1)), L' \ (L \ b(:, 2))];
uv(B, :) = Xb;

end

