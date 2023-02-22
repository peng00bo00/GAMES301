function [uv, M] = mean_value_nature_2(V, F, B, I, v1, v2)

nb = numel(B);
nv = size(V, 1);

[Laplacian] = mean_value(V, F);

%% boundary condition
% derivate of area at boundary
inv_len_b = reshape(1.0 ./ vecnorm(V(B, :) - V(B([2:nb 1]), :), 2, 2), [1 nb]);
% Mb = sparse([B B], nv + [B([nb 1:nb-1]) B([2:nb 1])], [ones(1, nb) -ones(1, nb)], nv * 2, nv * 2);
Mb1 = sparse([B B B], nv + [ B([nb 1:nb-1]) B([2:nb 1]) B],  [ inv_len_b([nb 1:nb-1]), -inv_len_b, inv_len_b-inv_len_b([nb 1:nb-1]) ], nv * 2, nv * 2);
Mb2 = sparse(nv + [B B B], [ B([nb 1:nb-1]) B([2:nb 1]) B], -[ inv_len_b([nb 1:nb-1]), -inv_len_b, inv_len_b-inv_len_b([nb 1:nb-1]) ], nv * 2, nv * 2);

M = Mb1 + Mb2;
M(1:nv, 1:nv) = Laplacian;
M((nv+1) : (nv*2), (nv+1) : (nv*2)) = Laplacian;

% fixed two points

% % find max distance
% [max_dis, ind] = max((V(B, 1) - V(B, 1).').^2 + (V(B, 2) - V(B, 2).').^2 + (V(B, 3) - V(B, 3).').^2, [], 'all', 'linear');
% [row, col] = ind2sub([numel(B), numel(B)], ind(1));
% % fprintf("err_dis = %.6f\n", max_dis - vecnorm(V(B(row), :) - V(B(col), :), 2).^2);
% v1 = B(row);
% v2 = B(col);

% v1 = B(1);
% v2 = B(idivide(int32(nb), 2) + 1);
% v2 = B(2);

Xb = zeros(nv, 2);
Xb(v1, :) = [1, 0];
Xb(v2, :) = [-1, 0];

%% solve

M_add = [M; sparse([1 2 3 4], [v1 v1 + nv v2 v2 + nv], 1, 4, nv * 2)];
y = zeros(nv * 2 + 4, 1);
y(nv * 2 + 1 : nv * 2 + 2) = Xb(v1, :);
y(nv * 2 + 3 : nv * 2 + 4) = Xb(v2, :);

uv_vec = M_add \ y;

% err = vecnorm(M_add * uv_vec - y, 'inf');
% fprintf("err_solver: %.10f\n", err);

uv = reshape(uv_vec, [nv, 2]);
end



