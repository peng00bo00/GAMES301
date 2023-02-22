function [I, B, Fin, Fsub] = local_patch(V, F, x, num)

nv = size(V, 1);

Gvv = sparse(F(:, [2 3 1]), F, 1, nv, nv);

vis = logical(zeros(1, nv));
que = zeros(1, num);

vis(x) = true;
que(1) = x;
s = 1;
t = 1;

while num > 0 && s <= t
    x = que(s);
    s = s + 1;
    
    ys = find(Gvv(:, x));
    
    for i=1:numel(ys)
        y = ys(i);
        if t < num && ~vis(y)
            vis(y) = true;
            que(t + 1) = y;
            t = t + 1;
        end
    end
end

Fin = vis(F(:, 1)) & vis(F(:, 2)) & vis(F(:, 3));
Fsub = F(Fin, :);
sub = find(vis);

[B, ~] = findBoundary(V, Fsub);
I = setdiff(sub, B);

end

