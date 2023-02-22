filename = 'mesh/cathead.obj';
file_big = 'mesh/camelhead.obj';
file_2d = 'mesh/2D.obj';
file_hilbert = 'mesh/ex.obj';
file_develop = 'mesh/CrumpledDevelopable.obj';
cow = 'mesh/cow.obj';

[v, f, ~, ~] = readObj(file_hilbert, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);

tic;
uv = conformal_nature(v, f, B, setdiff(1:nv, B), B(1), B(2));

% uv_lscm = conformal_nature(v, f, B, setdiff(1:nv, B));
% uv = mean_value_dirichlet(v, f, B, uv_lscm(B, :));

% uv = mean_value_nature(v, f, B, setdiff(1:nv, B));

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


