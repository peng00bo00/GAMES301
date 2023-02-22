function [V_cut, F_cut] = meshCut(V, F, w, b)
    nv = size(V, 1);

    in = (V * w - b) >= 0;
    f_have_v = in(F);
    f_in = f_have_v(:, 1) & f_have_v(:, 2) & f_have_v(:, 3);

    perm = zeros(nv, 1);
    inv_perm = find(in);
    perm(inv_perm) = (1:numel(inv_perm)).';
    
    V_cut = V(in, :);
    F_cut = perm(F(f_in, :));
end

