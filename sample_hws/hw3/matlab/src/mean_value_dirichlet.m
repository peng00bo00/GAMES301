function [uv] = mean_value_dirichlet(V, F, B, Bx)

nv = size(V, 1);
M = mean_value(V, F);
% M = cotLaplacian(V, F);

Xb = zeros(nv, 2);
Xb(B, :) = Bx;

y = -M * Xb;

I = setdiff(1:nv, B);
uv = zeros(nv, 2);
uv(I, :) = M(I, I) \ y(I, :);
uv = uv + Xb;

end

