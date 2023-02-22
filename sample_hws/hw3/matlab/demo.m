filename = 'mesh/cathead.obj';
file_big = 'mesh/camelhead.obj';
file_2d = 'mesh/2D.obj';
file_hilbert = 'mesh/ex.obj';
file_develop = 'mesh/CrumpledDevelopable.obj';
cow = 'mesh/cow.obj';

[v, f, ~, ~] = readObj(cow, false);
nv = size(v, 1);
[B, H] = findBoundary(v, f);
I = setdiff(1:nv, B);

[v1, v2] = find_diameter(v, B);
% v1 = B(1);
% v2 = B(206);

uv0 = conformal_nature(v, f, B, I, v1, v2);

flip_id = check_flip(v, f, uv0);
if (numel(flip_id) > 0)
    fprintf("\n\ninit flip: ");
    for i=1:numel(flip_id)
        fprintf("%d ", flip_id(i));
    end
    fprintf("\n");
end

uv1 = mean_value_dirichlet(v, f, B, uv0(B, :));

flip_id = check_flip(v, f, uv1);
if (numel(flip_id) > 0)
    fprintf("\n\nmean_value flip: ");
    for i=1:numel(flip_id)
        fprintf("%d ", flip_id(i));
    end
    fprintf("\n");
end

figure;
drawmesh(f, uv0, B);

figure;
drawmesh(f, uv1, B);


