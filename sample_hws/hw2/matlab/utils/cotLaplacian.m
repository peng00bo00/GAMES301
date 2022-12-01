function [L, Gcot] = cotLaplacian(V, F)

nv = size(V, 1);

% E = reshape((V(F(:, [2, 3, 1]), :) - V(F, :)).', [3 size(F)]);
% cots = squeeze(-dot(E(:, :, [2, 3, 1]), E(:, :, [3, 1, 2]), 1) ./ vecnorm(cross(E(:, :, [2, 3, 1]), E(:, :, [3, 1, 2]), 1), 2, 1));

E = reshape(V(F(:, [2, 3, 1]), :) - V(F, :), [size(F), 3]);
cots = -dot(E(:, [2, 3, 1], :), E(:, [3, 1, 2], :), 3) ./ vecnorm(cross(E(:, [2, 3, 1], :), E(:, [3, 1, 2], :), 3), 2, 3);

Gcot = sparse(F, F(:, [2, 3, 1]), cots, nv, nv);
L = Gcot + Gcot.';
L = diag(sparse(sum(L, 1))) - L;

end

% function [L] = cotLaplacian(X, T)
%     try 
%         L = evalin("base", cot_Laplacian);
%     catch
%         nv = size(X, 1);
%         Gvvn = sparse(T, T(:, [2, 3, 1]), T(:, [3, 1, 2]), nv, nv);
%         [ii, jj, kk] = find(Gvvn);
%         edge_cot = arrayfun(@(i, j, k) dot(X(i,:) - X(k,:), X(j,:) - X(k,:)) / norm(cross(X(i,:) - X(k,:), X(j,:) - X(k,:))), ii, jj, kk);
%         
%         L = sparse(ii, jj, edge_cot, nv, nv);
%         L = L + L.';
%     
%         L = diag(sum(L, 1)) - L;
%         assignin("base", 'cot_Laplacian', L);
%     end
% end

