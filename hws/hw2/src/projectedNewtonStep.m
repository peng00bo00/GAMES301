function [uv_new, stop] = projectedNewtonStep(V, F, uv, X1, As, PfPxs, maxIter, lam, c, tau)
%% One step projected Newton update
%% Args:
%%      V[nV, 3]: vertex coordinates in 3D
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%%      X1[3, 2, nF]: rest pose of each triangle
%%      As[nF, 1]: triangle areas
%%      PfPxs[4, 6, nF]: PfPx on different triangle
%%      maxIter: maximum iteration in line search
%%      lam: positive-definite parameter
%%      c: search control parameter
%%      tau: shrinkage factor in line search
%% Returns:
%%      uv_new[nV, 2]: updated uv coordinates
%%      stop: stop iteration flag

nV = size(V, 1);
nF = size(F, 1);

stop = false;

%% gradient
b = computeGradient(V, F, uv, X1, As, PfPxs);
g = reshape(b, [nV, 2]);

%% early stop
if norm(b) < 1e-4
    uv_new = uv;
    stop = true;
    return
end

%% Hessian
H = projectHessian(V, F, uv, X1, As, PfPxs);
H = H + lam * speye(nV * 2);

%% solve the update direction
d = -H \ b;
d = reshape(d, [nV, 2]);

%% one step update
[uv_new, E, alpha] = lineSearch(F, uv, X1, As, d, g, maxIter, c, tau);

fprintf("SD Energy: %.2e \tgradient: %.2e \talpha: %.2e\n", E, norm(b), alpha);

end