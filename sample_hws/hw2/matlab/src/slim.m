function [Xn, stop, M] = slim(X, mesh, e, wolfe_c1, wolfe_c2, lam, eps, chol_solver)

D = mesh.grad_operator;
sq_area = mesh.sq_area;
order = mesh.order;
nv = size(mesh.V, 1);
nf = size(mesh.F, 1);

grad_max = vecnorm(e.compute_gradient(X), 'inf');
% fprintf('grad_max = %.6f\n', grad_max);
if (grad_max <= eps) 
    Xn = X;
    stop = true;
    M = [];
    return;
else
    stop = false;
end

%%% Local Step   
Jvec = D * X;
Jblk = vec2block(Jvec);

[Wblk, Rblk] = local_step(Jblk);

Rvec = block2vec(Rblk);
Wvec = block2vec(Wblk);

%%% Global Step
tripI = [1:2*nf, 1:2*nf].';
tripI = [tripI; tripI + 2 * nf];
tripJ = [1:nf, 1:nf, nf+1:2*nf, nf+1:2*nf].';
tripJ = [tripJ; tripJ + 2 * nf];

W = sparse(tripI, tripJ, [Wvec; Wvec], nf * 4, nf * 4);

A = sq_area * W * D;
b = sq_area * (W * Rvec);

% unit test: check residual between grad_slim and grad_sd
% Mo = A.' * A;
% yo = A.' * b;
% grad_slim = full(Mo * X - yo);
% grad_real = e.compute_gradient(X);
% 
% residual = vecnorm((grad_slim - grad_real), 'inf');
% fprintf('residual between grad_slim and grad_real: %.6f\n', residual);

M = A.' * A + speye(nv * 2) .* lam;
y = A.' * b + X .* lam;

p = zeros(nv * 2, 1);

if chol_solver
    L = chol(M(order, order), 'lower');
    p(order) = (L.' \ (L \ y(order)));
else
    p(order) = M(order, order) \ y(order);
end

dx = (p - X); % descent direction

%%% Line Search
min_step_to_singularity = compute_min_step_to_singularity(D, X, dx);
[Xn, alpha] = line_search(X, dx, min_step_to_singularity, wolfe_c1, wolfe_c2, e, min([1.0, 0.8 * min_step_to_singularity]));

fprintf('2-norm of dx: %.6f,\talpha: %.6f,\tstep_to_singularity: %.6f\n', vecnorm(dx, 2), alpha, min_step_to_singularity);

end

function B = vec2block(A)
    % input:  A (nf * 4, 1) [J_11, J_21, J_12, J_22]
    % output: B (2, 2, nf)
    nf = size(A, 1) / 4;
    B = reshape(reshape(A, [nf 4]).', [2 2 nf]);
end

function A = block2vec(B)
    % input:  B (2, 2, nf)
    % output: A (nf * 4, 1) [J_11, J_21, J_12, J_22]
    nf = size(B, 3);
    A = reshape(reshape(B, [4 nf]).', [4 * nf, 1]);
end