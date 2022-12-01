addpath("io\");
addpath("mesh\");
addpath("src\");
addpath("utils\");

file_camel = 'mesh/camelhead.obj';

[v, f, ~, ~] = readObj(file_camel, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);
nf = size(f, 1);

%% parameter
eps = 1e-4;
lam = 1e-4;
wolfe_c1 = 5e-6;
wolfe_c2 = 0.9;

%% initialize
uv = conformal(v, f, B);
% uv = tutte(v, f, B);
X = reshape(uv, [nv * 2, 1]);

flip_id = check_flip(v, f, uv);
if (numel(flip_id) > 0)
    fprintf('\n\twarning: %d triangles flipped!!\n\n', numel(flip_id));
end

%% grad_operator
[grad_2d, ~, Area] = deform_grad_operator(v, f);
D1 = grad_2d(1:2:nf*2, :);
D2 = grad_2d(2:2:nf*2, :);

D = [D1 sparse(nf, nv); 
     sparse(nf, nv) D1;
     D2 sparse(nf, nv);
     sparse(nf, nv) D2];

%% sq_area
sq_area = diag(sparse(repmat(sqrt(Area), [4, 1])));

%% reordering
L = cotLaplacian(v, f);
p = dissect(repmat(L, [2, 2]));

%% energy
e = SymDirichlet(v, f);

%% store them into mesh struct
mesh.V = v;
mesh.F = f;
mesh.grad_operator = D;
mesh.sq_area = sq_area;
mesh.order = p;

figure;
h = drawmesh(f, uv, B);

fprintf("step 0, energy: %.4f\n", e.compute_energy(reshape(uv, [nv * 2, 1])));

for i = 1:100
    fprintf('\nstep %d\n', i);

    tic;
    [X, stop] = proj_newton(X, mesh, e, wolfe_c1, wolfe_c2, lam, eps, true);
%     [X, stop] = slim(X, mesh, e, wolfe_c1, wolfe_c2, lam, eps, true);
    toc;

    if stop
        fprintf("converge after step %d!!\n", i - 1);
        break;
    end

    %% eval energy
    energy = e.compute_energy(X);
    grad = vecnorm(e.compute_gradient(X), 'inf');
    fprintf("energy: %.6f,\tgradient: %.6f\n", energy, grad);

    uv = reshape(X, [nv, 2]);

    %% check flip
    flip_id = check_flip(v, f, uv);
    if (numel(flip_id) > 0)
        fprintf('\n\twarning: %d triangles flipped!!\n\n', numel(flip_id));
    end

    %% draw
    set(h(1), 'vertices', uv);
    set(h(2), 'XData', uv([B B(1)], 1), 'YData', uv([B B(1)], 2));
    drawnow;

    pause(0.01);
end