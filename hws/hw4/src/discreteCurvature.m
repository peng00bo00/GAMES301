function [K, kappa] = discreteCurvature(V, F, B)
%% Compute discrete Gaussian and geodesic curvature
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%%      B[1, nB]: boundary vertex index
%% Returns:
%%      K[nV, 1]: discrete Gaussian curvature
%%      kappa[nB, 1]: discrete geodesic curvature at the boundary

nV = size(V, 1);
nF = size(F, 1);

%% edges
Es    = reshape(V(F(:, [2, 3, 1]), :) - V(F, :), [size(F), 3]);
Enorm = vecnorm(Es, 2, 3);
Edir  = Es ./ Enorm;

%% interior angles
coss   =-dot(Edir(:, :, :), Edir(:, [3, 1, 2], :), 3);
angles = acos(coss);

%% sum over vertices
A = sparse(repmat((1:nF)', 1, 3), F, angles, nF, nV);
angles = full(sum(A, 1))';

%% discrete curvatures
K     = 2*pi - angles;
kappa = pi - angles(B);

%% set K=0 at the boundary
K(B) = 0;

end