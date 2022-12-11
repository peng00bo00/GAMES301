function [uv_new, Enew, alpha] = lineSearch(F, uv, X1, As, d, g, maxIter, c, gamma)
%% Line search to update parameterization

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

    alpha = alpha * gamma;
end

% c1 = 5e-6;
% c2 = 0.9;
% 
% nV = size(V, 1);
% 
% alpha_max = 0.99*initStepSize(F, uv, d);
% alpha_min = 0.0;
% alpha     = min([0.8*alpha_max, 1.0]);
% 
% Eold = computeTotalSDEnergy(F, uv, X1, As);
% Gold = sum(d .* g, "all");
% 
% uv_new = uv + alpha * d;
% Enew = computeTotalSDEnergy(F, uv_new, X1, As);
% 
% for i=1:maxIter
%     if Enew > Eold + c1*alpha*Gold
%         alpha_max = alpha;
%         alpha = (alpha_min + alpha_max) / 2;
%     else
%         grad = computeGradient(V, F, uv, X1, As, PfPxs);
%         grad = reshape(grad, [nV, 2]);
% 
%         Gnew = sum(d .* grad, "all");
% 
%         if abs(Gnew) > c2 * abs(Gold)
%             if Gnew > 0
%                 alpha_max = alpha;
%             else
%                 alpha_min = alpha;
%             end
% 
%             alpha = (alpha_min + alpha_max) / 2;
%         else
%             return
%         end
%     end
% 
%     uv_new = uv + alpha * d;
%     Enew = computeTotalSDEnergy(F, uv_new, X1, As);
% end

end

function alpha = initStepSize(F, uv, d)
%% A helper function to find initial step size in line search
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

M = D * matrixInv2x2(Ds);
[~, S, ~] = svd_rv(M);

t1 = max(diag(S));
if t1 <= 0
    return
else
    alpha = 1 / t1;
end

% %% quadratic equation
% detDs= det(Ds);
% detD = det(D);

% A = detD;
% B = detDs + detD - det(Ds-D);
% C = detDs;
% 
% delta = B*B - 4*A*C;
% 
% if delta < 0
%     return
% else
%     alpha = (-B-sqrt(delta))/(2*A);
%     if alpha < 0
%         alpha = inf;
%     end
% end

end