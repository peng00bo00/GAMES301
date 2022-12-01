function [Xn, stop, H] = proj_newton(X, mesh, e, wolfe_c1, wolfe_c2, lam, eps, chol_solver)

nv = size(mesh.V, 1);
nf = size(mesh.F, 1);

b = -e.compute_gradient(X);

grad_max = vecnorm(b, 'inf');
% fprintf('grad_max = %.6f\n', grad_max);
if (grad_max <= eps)
    stop = true;
    Xn = X;
    H = [];
    return;
else
    stop = false;
end

H = e.compute_hessian(X);
M = H + lam * speye(nv * 2);

p = mesh.order;

d = zeros(nv * 2, 1);
if chol_solver
    L = chol(M(p, p), 'lower');
    d(p) = L.' \ (L \ b(p));
else
    d(p) = M(p, p) \ b(p);
end

step_max = compute_min_step_to_singularity(mesh.grad_operator, X, d);
% alpha = min([1.0, 0.5 * step_max]);
% Xn = X + alpha * d;
Xn = line_search(X, d, step_max, wolfe_c1, wolfe_c2, e, min([1, 0.8 * step_max]));

% len = vecnorm(d, 2);
fprintf('len of direction = %.8f,\tstep_max = %.8f\n', vecnorm(d, 2), step_max);

end

