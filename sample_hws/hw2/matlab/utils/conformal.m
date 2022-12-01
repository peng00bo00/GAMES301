function [uv] = conformal(V, F, B)

nb = numel(B);
nv = size(V, 1);

M = cotLaplacian(V, F);
Gcot = -M + diag(sparse(diag(M)));

%% boundary condition
% v1 = B(1);
% v2 = B(mod(int32(nb), 2) + 1);
% I = setdiff(1:nv, [v1, v2]);
% Xb = zeros(nv, 2);
% Xb(v1, :) = [-1, 0];
% Xb(v2, :) = [1, 0];

I = setdiff(1:nv, B);
Xb = zeros(nv, 2);
Xb(B, :) = setBoundary(nb, 'circle');
%% solve

y = Gcot * sparse(Xb);

uv = zeros(nv, 2);
uv(I, :) = M(I, I) \ y(I, :);

% L = chol(M(I, I), 'lower'); 
% uv(I, :) = L.' \ (L \ y(I, :));

uv = uv + Xb;
end

