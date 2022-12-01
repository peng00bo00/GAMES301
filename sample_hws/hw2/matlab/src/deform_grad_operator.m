function [grad_2d, grad_3d, Area] = deform_grad_operator(V, F)

nv = size(V, 1);
nf = size(F, 1);

Edge = reshape(V(F(:, [2 3 1]), :) - V(F, :), [nf 3 3]);
NormW = squeeze(cross(Edge(:, 1, :), Edge(:, 2, :), 3)); % [f, 3]
Area = vecnorm(NormW, 2, 2); % [f, 1]
Norm = NormW ./ Area;   % [f, 3]
Norm3 = repmat(reshape(Norm, [nf 1 3]), [1 3 1]); % [f, 3, 3]

Gb = cross(Norm3, Edge(:, [2 3 1], :), 3) ./ Area;
% tripI(i, j, k) = 3 * (i - 1) + k
tripI = repmat(reshape(reshape(1:3*nf, [3 nf]).', [nf 1 3]), [1 3 1]);
tripJ = repmat(F, [1 1 3]);

grad_3d = sparse(reshape(tripI, [nf*3*3,1]), reshape(tripJ, [nf*3*3,1]), reshape(Gb,[nf*3*3,1]), nf * 3, nv);

% compute local frame basis

lx = squeeze(Edge(:, 1, :));
lx = lx ./ vecnorm(lx, 2, 2);   % local frame base x
ly = cross(Norm, lx, 2);        % local frame base y

JJ = reshape(1:nf*3, [3 nf]).';
II = repmat((1:2:nf*2).', [1 3]);
projLocal = sparse([II; II+1], [JJ; JJ], [lx; ly], nf * 2, nf * 3);

grad_2d = projLocal * grad_3d;

end

