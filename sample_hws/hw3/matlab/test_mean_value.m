file_hilbert = 'mesh/ex.obj';

[v, f, ~, ~] = readObj(file_hilbert, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);

tic;
uv = mean_value_nature_2(v, f, B, setdiff(1:nv, B), B(1), B(2));
toc;

flip_id = check_flip(v, f, uv);
if (numel(flip_id) > 0)
    for i=1:numel(flip_id)
        fprintf("%d flipped!!\n", flip_id(i));
    end
end

figure;
drawmesh(f, v, B);

figure;
drawmesh(f, uv, B);


