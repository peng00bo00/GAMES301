filename = 'mesh/cathead.obj';
file_big = 'mesh/camelhead.obj';
file_2d = 'mesh/2D.obj';
file_hilbert = 'mesh/ex.obj';
file_develop = 'mesh/CrumpledDevelopable.obj';
bunny_head = 'mesh/Bunny_head.obj';
tri = 'mesh/square_sparse.obj';
fertility = 'mesh/Fertility.obj';
vest = 'mesh/vest.obj';

[v, f, ~, ~] = readObj(bunny_head, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);
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

for k=1:5
    vv = [uv0, zeros(nv, 1)];
    uv1 = mean_value_nature(vv, f, B, I, v1, v2);

    flip_id = check_flip(v, f, uv1);
    if (numel(flip_id) > 0)
        fprintf("\n\ncurrent flip: ");
        for i=1:numel(flip_id)
            fprintf("%d ", flip_id(i));
        end
        fprintf("\n");
    else
        fprintf("\nno flip!!\n");
    end

    uv0 = uv1;
end

figure;
drawmesh(f, uv0, B);


