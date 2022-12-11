function [uv_new, stop] = projectedNewtonStep(V, F, uv, X1, As, PfPxs, maxIter, lam, c, gamma)
%% One step projected Newton solver
%% Args:
%%      V[nV, 3]: vertex coordinates in 3D
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%%      maxIter: maximum iteration in line search
%%      lam: positive-definite parameter
%%      c: search control parameter
%%      gamma: shrinkage factor in line search
%% Returns:
%%      uv_new[nV, 2]: updated uv coordinates
%%      stop: stop iteration flag

nV = size(V, 1);
nF = size(F, 1);

stop = false;

%% Newton's method
b = computeGradient(V, F, uv, X1, As, PfPxs);
g = reshape(b, [nV, 2]);

%% early stop
if norm(b) < 1e-4
    uv_new = uv;
    stop = true;
    return
end

H = projectHessian(V, F, uv, X1, As, PfPxs);
H = H + lam * speye(nV * 2);

%% solve update direction
d = -H \ b;
d = reshape(d, [nV, 2]);

%% one step update
[uv_new, E, alpha] = lineSearch(F, uv, X1, As, d, g, maxIter, c, gamma);

fprintf("SD Energy: %.2e \tgradient: %.2e \talpha: %.2e\n", E, norm(b), alpha);

end