function [uv] = conformal_nature_multi(V, F, B, H, roi, v1, v2)

nb = numel(B);
nv = size(V, 1);

Laplacian = cotLaplacian(V, F);

%% boundary condition
% derivate of area at boundary
Mb = sparse([B B], nv + [B([nb 1:nb-1]) B([2:nb 1])], [ones(1, nb) -ones(1, nb)], nv * 2, nv * 2);
M = Mb + Mb.';

for i=1:numel(H)
    cur_B = H{i};
    cur_nb = numel(cur_B);
    Mb = sparse([cur_B cur_B], nv + [cur_B([cur_nb 1:cur_nb-1]) cur_B([2:cur_nb 1])], [ones(1, cur_nb), -ones(1, cur_nb)], nv * 2, nv * 2);
    M = M + Mb + Mb.';
end

M(1:nv, 1:nv) = Laplacian;
M((nv+1) : (nv*2), (nv+1) : (nv*2)) = Laplacian;

% fixed two points

roi = setdiff([roi, roi + nv], [v1, v2, v1 + nv, v2 + nv]);

Xb = zeros(nv, 2);
Xb(v1, :) = [1, 0];
Xb(v2, :) = [-1, 0];

%% solve

uv_vec = zeros(nv * 2, 1);

y = -M * sparse(reshape(Xb, [nv * 2, 1]));

% L = chol(M(roi, roi), 'lower'); 
% uv_vec(roi) = L.' \ (L \ y(roi));

uv_vec(roi) = M(roi, roi) \ y(roi);
uv_vec = uv_vec + reshape(Xb, [nv * 2, 1]);

uv = reshape(uv_vec, [nv, 2]);
end


