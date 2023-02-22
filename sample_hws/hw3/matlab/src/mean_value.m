function [W, len] = mean_value(V, F)

nv = size(V, 1);

E = reshape(V(F(:, [2 3 1]), :) - V(F, :), [size(F), 3]);
% E = reshape(V(F(:, [2 3 1]), :) - V(F, :), [size(F), 3]);
len = vecnorm(E, 2, 3);
cos = -dot(E, E(:, [3 1 2], :), 3) ./ (len .* len(:, [3 1 2]));

tan_hf = sqrt((1-cos)./(1+cos));
Wpos = sparse(F, F(:, [2 3 1]), tan_hf ./ len, nv, nv);
Wneg = sparse(F(:, [2 3 1]), F, tan_hf(:, [2 3 1]) ./ len, nv, nv);
W = Wpos + Wneg;

if (numel(find(W<0, 1)) > 0)
    fprintf("warning: mean_value_coordinates negative!!\n");
end

W = diag(sparse(sum(W, 2))) - W;

end

