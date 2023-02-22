function [v1, v2] = find_diameter(V, B)

% find max distance
[max_dis, ind] = max((V(B, 1) - V(B, 1).').^2 + (V(B, 2) - V(B, 2).').^2 + (V(B, 3) - V(B, 3).').^2, [], 'all', 'linear');
[row, col] = ind2sub([numel(B), numel(B)], ind(1));
% fprintf("err_dis = %.6f\n", max_dis - vecnorm(V(B(row), :) - V(B(col), :), 2).^2);
v1 = B(row);
v2 = B(col);

end

