function [eigs, lams] = SDEigenSystem(F)
%% Analytic eigensystem of symmetric Dirichlet energy
%% Args:
%%      F[2, 2]: deformation gradient
%% Returns:
%%      eigs[4, 1]: eigen values
%%      lams[4, 4]: each column is an eigen vector

%% SVD
[R, S, U, Sigma, V] = polarSVD(F);

I1 = trace(S);
I2 = norm(S, "fro")^2;
I3 = det(S);

%% eigen values
eigs = zeros(4, 1);
eigs(1) = 1 + 3/Sigma(1,1)^4;
eigs(2) = 1 + 3/Sigma(2,2)^4;
eigs(3) = 1 + 1/I3^2 + I2/I3^3;
eigs(4) = 1 + 1/I3^2 - I2/I3^3;

%% eigen vectors
lams = zeros(4);

D1 = U * [1, 0; 0, 0] * V';
D2 = U * [0, 0; 0, 1] * V';
L  = U * [0, 1; 1, 0] * V' / sqrt(2);
T  = U * [0,-1; 1, 0] * V' / sqrt(2);

lams(:, 1) = vec(D1);
lams(:, 2) = vec(D2);
lams(:, 3) = vec(L);
lams(:, 4) = vec(T);

end