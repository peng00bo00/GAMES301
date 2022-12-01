function [flip_id] = check_flip(V, F, uv)

    nf = size(F, 1);

    D = deform_grad_operator(V, F);
    J = full(D * uv);
    J = reshape(J.', [2 2 nf]);

    sgn = J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :);
    reshape(sgn, [1 nf]);
    
    flip_id = find(sgn(sgn < 0));

end

