clear; clc;

addpath("io\");
addpath("mesh\");
addpath("src\");
addpath("utils\");

filename = 'mesh/cathead.obj';
file_camel = 'mesh/camelhead.obj';
file_face = 'mesh/Nefertiti_face.obj';
file_jury = 'mesh/camelhead-proj.obj';
file_cow = 'mesh/cow.obj';
file_ex = 'mesh/ex.obj';

[v, f, ~, ~] = readObj(file_ex, false);
% [vv, ff, vt, ~] = readObj(file_jury, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);
nf = size(f, 1);

%% parameter
eps = 1e-8;
lam = 1e-4;
wolfe_c1 = 5e-6;
wolfe_c2 = 0.9;

%% initialize
uv = conformal(v, f, B);
% % uv = tutte(v, f, B);
% 
% flip_id = check_flip(v, f, uv);
% if numel(flip_id) > 0
%     fprintf('\n\twarning: initial %d triangle flipped!!\n\n', numel(flip_id));
% end
% 
% X = reshape(uv, [nv * 2, 1]);
% 
% %% grad_operator
% [grad_2d, ~, Area] = deform_grad_operator(v, f);
% D1 = grad_2d(1:2:nf*2, :);
% D2 = grad_2d(2:2:nf*2, :);
% 
% D = [D1 sparse(nf, nv); 
%      sparse(nf, nv) D1;
%      D2 sparse(nf, nv);
%      sparse(nf, nv) D2];
% 
% %% sq_area
% sq_area = diag(sparse(repmat(sqrt(Area), [4, 1])));
% 
% %% reordering
% L = cotLaplacian(v, f);
% p = dissect(repmat(L, [2, 2]));
% 
% %% energy
% e = SymDirichlet(v, f);
% 
% %% store them into mesh struct
% mesh.V = v;
% mesh.F = f;
% mesh.grad_operator = D;
% mesh.sq_area = sq_area;
% mesh.order = p;
% 
% figure;
% h = drawmesh(f, uv, B);
% 
% fprintf("step 0, energy: %.4f\n", e.compute_energy(reshape(uv, [nv * 2, 1])));
% 
% XX = zeros(nv * 2, 1);
% 
% 
% for i = 1:1000
%     fprintf('\nstep %d\n', i);
% 
%     if i == 1
%         [X, stop] = proj_newton(X, mesh, e, wolfe_c1, wolfe_c2, 1e-4, eps, false);
%     elseif i <= 25
%         [X, stop] = slim(X, mesh, e, wolfe_c1, wolfe_c2, 1e-5, eps, false);
%     else
%         [X, stop] = proj_newton(X, mesh, e, wolfe_c1, wolfe_c2, 1e-8, eps, true);
%     end
% 
%     if stop
%         break;
%     end
% 
%     %% show
%     energy = e.compute_energy(X);
%     grad = vecnorm(e.compute_gradient(X), 'inf');
%     fprintf("energy: %.6f,\tgradient: %.6f\n", energy, grad);
% 
%     uv = reshape(X, [nv, 2]);
% 
%     % check flip
%     flip_id = check_flip(v, f, uv);
%     if (numel(flip_id) > 0)
%         fprintf('\n\twarning: %d triangles flipped!!\n\n', numel(flip_id));
%     end
% 
% 
%     set(h(1), 'vertices', uv);
%     set(h(2), 'XData', uv([B B(1)], 1), 'YData', uv([B B(1)], 2));
%     drawnow;
% 
%     pause(0.01);
% end
% 
% 
% % std = norm(vt, "inf")
% % err = norm(uv - vt, "inf")
% % 
% % drawmesh(f, uv, B);
% % 
% % e_check = SymDirichlet([vt, zeros(nv, 1)], f);
% % err_energy = e_check.compute_energy(reshape(uv, [nv * 2, 1]));
% % 
% % fprintf('energy to std; %.4f\n', err_energy);
% 
% 
% %     if i == 5
% %         [X, stop] = proj_newton(X, mesh, e, wolfe_c1, wolfe_c2, 1e-6, false);
% %     elseif i <= 35
% %         [X, stop] = slim(X, mesh, e, wolfe_c1, wolfe_c2, 1e-3, false);
% %     elseif i <= 80
% %         [X, stop] = proj_newton(X, mesh, e, 0.15, 1.1, 1e-7, false);
% %     elseif mod(i, 12) == 0
% %         [X, stop] = proj_newton(X, mesh, e, -0.12, 1.75, 1e-8, false);
% %     else
% %         [X, stop] = proj_newton(X, mesh, e, 0.005, 1.75, 1e-8, false);
% %     end
% 
% %     if i <= 1
% %         X = proj_newton(X, mesh, e, 0.1, 0.99, 1e-5, false);
% %     elseif i <= 18
% %         X = slim(X, mesh, e, 0.01, 0.99, 1e-6, false);
% %     else
% %         X = proj_newton(1.02 * X - 0.02 * XX, mesh, e, -0.5, 1.75, 1e-8, true);
% %     end
% % 
% %     XX = X;