filename = 'mesh/cathead.obj';
file_big = 'mesh/camelhead.obj';
file_2d = 'mesh/2D.obj';
file_hilbert = 'mesh/ex.obj';
file_develop = 'mesh/CrumpledDevelopable.obj';
bunny_head = 'mesh/Bunny_head.obj';
tri = 'mesh/square_sparse.obj';
fertility = 'mesh/Fertility.obj';
vest = 'mesh/vest.obj';

[v, f, ~, ~] = readObj(file_develop, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);

[Fsub, I1, B1, H1] = make_holes(v, f, 9, 40);
roi = [I1, B1];
for i=1:numel(H1)
    roi = [roi, H1{i}];
end

I2 = I1;
for i=1:numel(H1)
    I2 = [I2, H1{i}];
end

figure;
drawmesh(Fsub, v, B);

% text(v(:, 1), v(:, 2), v(:, 3), num2cell(1:nv));
% plot3(v(roi, 1), v(roi, 2), v(roi, 3), 'ro');

nv = size(v, 1);

v1 = 5;
v2 = 239;

uv0 = conformal_nature_multi(v, Fsub, B1, H1, roi, v1, v2);
% uv0 = conformal_nature(v, f, B, setdiff(1:nv, B), v1, v2);

flip_id = check_flip(v, Fsub, uv0);
if (numel(flip_id) > 0)
    fprintf("\n\ninit flip: ");
    for i=1:numel(flip_id)
        fprintf("%d ", flip_id(i));
    end
    fprintf("\n");
end

figure;
drawmesh(Fsub, uv0, B);


