function uv = MVCTutte(V, F)
%% Tutte's embedding with mean value coordinates
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      uv[nV, 2]: uv coordinates

nV = size(V, 1);

[B, ~] = findBoundary(V, F);
nB = length(B);

%% LSCM to pin the boundary
uv = LSCM(V, F);

%% edges
E     = reshape(V(F(:, [2, 3, 1]), :) - V(F, :), [size(F), 3]);
Enorm = vecnorm(E, 2, 3);
Edir  = E ./ Enorm;

%% trigonometry
coss =-dot(Edir(:, :, :), Edir(:, [3, 1, 2], :), 3);
sins = vecnorm(cross(Edir(:, :, :), -Edir(:, [3, 1, 2], :), 3), 2, 3);

%% tan(theta/2) = sin(theta) / (1+cos(theta))
tanHalf = sins ./ (coss + 1);

%% set up linear system
tripI = [F F];
tripJ = [F(:, [2 3 1]) F(:, [3 1 2])];
tripV = [tanHalf ./ Enorm(:,:) tanHalf ./ Enorm(:, [3 1 2])];

A = sparse(tripI, tripJ, tripV, nV, nV);
A = A - diag(sum(A, 2));

%% pinned boundary points
A(B, :) = sparse(nB, nV);
A(B, B) = speye(nB);

b = zeros(nV, 2);
b(B,:) = uv(B, :);

%% solve linear system
uv = A \ b;

end