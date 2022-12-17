function uv = LSCM(V, F)
%% LSCM parameterization
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      uv[nV, 2]: uv coordinates

nV = size(V, 1);
nF = size(F, 1);

%% areas
AT = doubleArea(V, F);
AT_sq = sqrt(AT);

%% rest pose
Xs = zeros(nF, 3, 2);
for i=1:nF
    Xs(i,:,:) = project2Plane(V(F(i, :), :));
end

%% set up sparse complex matrix
A = (Xs(:, [3 1 2], 1) - Xs(:, [2 3 1], 1)) ./ AT_sq;   %% real part
B = (Xs(:, [3 1 2], 2) - Xs(:, [2 3 1], 2)) ./ AT_sq;   %% imag part

MA = sparse(repmat((1:nF)', 1, 3), F, A, nF, nV);
MB = sparse(repmat((1:nF)', 1, 3), F, B, nF, nV);

%% pin 2 points on the boundary
[~, ~, T] = pinBoundary(V, F);

%% swap columns, now the pinned points are the last 2 ones
MA = MA * T; MB = MB * T;

%% split free and pinned vertices
Af = MA(:, 1:nV-2); Ap = MA(:, end-1:end);  %% real part
Bf = MB(:, 1:nV-2); Bp = MB(:, end-1:end);  %% imag part

AM = [Af -Bf; Bf Af];
b  =-[Ap -Bp; Bp Ap] * [0; 1; 0; 0];

%% solve linear system
uv = AM \ b;
uv = reshape(uv, [nV-2 2]);

%% fix pinned points
uv = [uv; [0 0]; [1 0]];

%% transform back to original vertex index
uv = T * uv;

end

function [b1, b2, T] = pinBoundary(V, F)
%% A helper function to pin 2 points on the boundary
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      b1: pinned boundary vertex
%%      b2: pinned boundary vertex
%%      T[nV, nV]: vertex index transform matrix

nV = size(V, 1);

[B, ~] = findBoundary(V, F);
nB = length(B);

% %% distance between vertices
% VB = V(B, :);
% XY = VB * VB'; X2 = diag(XY);
% D = X2 + X2' - 2*XY;
% D = triu(D);
% 
% %% select the farthest pair
% [i j] = find(ismember(D, max(D(:))));
% b1 = B(i); b2 = B(j);

%% select the first and middle boundary points
b1 = B(1); b2 = B(round(nB/2));

%% sparse index transform matrix
tripI = 1:nV; tripJ = 1:nV; tripV = ones(1, nV);
tripJ([b1 end-1]) = tripJ([end-1 b1]);
tripJ([b2 end])   = tripJ([end b2]);

T = sparse(tripI, tripJ, tripV, nV, nV);

end