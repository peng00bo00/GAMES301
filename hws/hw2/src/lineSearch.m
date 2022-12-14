function [uv_new, Enew, alpha] = lineSearch(F, uv, X1, As, d, g, maxIter, c, tau)
%% Line search to update parameterization
%% Args:
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%%      X1[3, 2, nF]: rest pose of each triangle
%%      As[nF, 1]: triangle areas
%%      d[nV, 2]: update direction
%%      g[nV, 2]: gradients
%%      maxIter: maximum iteration
%%      c: search control parameter
%%      tau: shrinkage factor in line search
%% Returns:
%%      uv_new[nV, 2]: updated uv coordiantes
%%      Enew: updated energy
%%      alpha: optimal step size

%% initialize step size with flip check
alpha = initStepSize(F, uv, d);
alpha = min(0.99*alpha, 1.0);

%% backtracking line search
Eold = computeTotalSDEnergy(F, uv, X1, As);
Enew = Eold;

uv_new = uv;

for i=1:maxIter
    uv_new = uv + alpha * d;
    Enew = computeTotalSDEnergy(F, uv_new, X1, As);

    %% Armijo–Goldstein condition
    if Enew <= Eold + c * alpha * sum(d .* g, "all")
        return;
    end

    alpha = alpha * tau;
end

end


function alpha = initStepSize(F, uv, d)
%% A helper function to find initial step size in line search
%% Args:
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%%      d[nV, 2]: update direction
%% Returns:
%%      alpha: maximum allowed step size

alpha = inf;
nF = size(F, 1);

for i=1:nF
    x = uv(F(i, :), :);
    dx= d(F(i, :), :);
    
    alpha = min([alpha, flipCheckStepSize(x, dx)]);
end

end


function alpha = flipCheckStepSize(x, dx)
%% A helper function to find maximum step size to prevent flip
%% Args:
%%      x[3, 2]: triangle coordinates in 2D, each row is a coordinate vector
%%      dx[3, 2]: update direction
%% Return:
%%      alpha: maximum step size

alpha = inf;

Ds= x(2:end, :) - x(1, :);
D = dx(2:end, :) - dx(1, :);

%% quadratic equation
detDs= det(Ds);
detD = det(D);

A = detD;
B = detDs + detD - det(Ds-D);
C = detDs;

delta = B*B - 4*A*C;

if delta < 0
    return
else
    alpha = (-B-sqrt(delta))/(2*A);
    if alpha < 0
        alpha = inf;
    end
end

end