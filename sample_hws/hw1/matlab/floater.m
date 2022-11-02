function [uv] = floater(X, F, B)

% input: X (geometry), F (topology), B (boundary indices)
% output: parameterization for each vertex

nv = size(X, 1);
nb = numel(B);
nf = size(F, 1);

% calculate edge length and angle of each triangle
edge = reshape(X(F(:, [2, 3, 1]), :) - X(F(:, :), :), [nf, 3, 3]);
len_t = sqrt(sum(edge.^2, 3));
angle_t = acos(-dot(edge(:, :, :), edge(:, [3, 1, 2], :), 3) ./ (len_t.*len_t(:, [3, 1, 2])));

% Gvvn(j, i) = k, where k is another vertex on triangle including edge (i, j)
% Gvvi(j, i) is an index of edge (i, j) on F
Gvvn = sparse(F(:, [2, 3, 1]), F, F(:, [3, 1, 2]), nv, nv);
Gvvi = sparse(F(:, [2, 3, 1]), F, 1:numel(F), nv, nv);

D = full(sum(spones(Gvvn), 1));  % degree of each vertex
Interior = setdiff(1:nv, B);

function [J, ind] = sorting(i, deg)
    J = zeros(deg, 1);
    Gcol = Gvvn(:, i);
    J(1) = find(Gcol, 1);
    for k=2:deg
        J(k) = Gcol(J(k - 1));
    end
    % ind denotes indices of elements on 1-ring, used to get length and angle on 1-ring later. 
    ind = full(Gvvi(J, i));
end

count = D(Interior);

[Jcell, indcell] = arrayfun(@sorting, Interior, count, 'UniformOutput', false);

inds = vertcat(indcell{:});
angles = angle_t(inds); % angle on 1-ring
lens = len_t(inds);     % length on 1-ring

% we use count to blocking(mat2cell) angles and lens to cell, each cell contains
% information for calculating floater weights on 1-ring.
ijvcell = cellfun(@calc_floater_weights, num2cell(Interior).', Jcell.', mat2cell(angles, count), mat2cell(lens, count), 'UniformOutput', false);

ijvs = vertcat(ijvcell{:});
W_floater = sparse(ijvs(:, 1), ijvs(:, 2), ijvs(:, 3), nv, nv); % sparse matrix for floater weights.

Xb = zeros(nv, 2);
Xb(B, :) = setBoundary(nb, 'circle');

y = W_floater * Xb;
M = speye(nv - nb) - W_floater(Interior, Interior);

uv = zeros(nv, 2);
uv(Interior, :) = M \ y(Interior, :);
uv = uv + Xb;

end

function [ijv] = calc_floater_weights(i, id, ang, len)

% i : vertex id
% id : global indices of elements on 1-ring
% ang: angle of each element on 1-ring
% len: length of each out-going edge on 1-ring

% ijv: triplets for constructing sparse floater matrix

    deg = numel(id);

    ang = (ang ./ sum(ang)) * 2 * pi;   % normalize angle
    ang = ang([deg 1:deg-1]);           % shift
    ang = cumsum(ang);                  % prefix sum

    % Turn 720 degree along the 1-ring, calculate the prefix sum of the angles, and we get double_ang
    double_ang = repmat(ang, [2, 1]);
    double_ang(deg+1:2*deg) = double_ang(deg+1:2*deg) + 2 * pi;

    % By calculating difference between the prefix sum of angle, we get the interval sum. 
    % For each vertex i on 1-ring, we find the first vertex j that turns
    % more than 180 degrees relative to i, then we find the triangle
    % containing the (0,0).
    [J, I] = find(diff((double_ang - ang.') >= pi));

    % tri_id contains local indices of triangle containing (0,0)
    tri_id = [I, J, J + 1];
    tri_id(tri_id >= deg+0.5) = tri_id(tri_id >= deg+0.5) - deg; % mod it to [1, deg]

    % area = 1/2 * a * b * sin\theta, using area to calculate Barycentric
    % coordinates。
    area = sin(ang(tri_id(:, [3, 1, 2])) - ang(tri_id(:, [2, 3, 1]))) .* len(tri_id(:, [3, 1, 2])) .* len(tri_id(:, [2, 3, 1]));
    area = area ./ sum(area, 2);
    area = area ./ deg;

    num = numel(tri_id);
    ijv = [ones(num, 1) .* i, id(reshape(tri_id, [num, 1])), reshape(area, [num, 1])];
end
