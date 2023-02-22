filename = 'mesh/cathead.obj';
file_big = 'mesh/camelhead.obj';
file_2d = 'mesh/2D.obj';
file_hilbert = 'mesh/ex.obj';
file_develop = 'mesh/CrumpledDevelopable.obj';
bunny_head = 'mesh/Bunny_head.obj';
tri = 'mesh/square_sparse.obj';
fertility = 'mesh/Fertility.obj';
vest = 'mesh/vest.obj';

[v, f, ~, ~] = readObj(vest, false);
[B, H] = findBoundary(v, f);

nv = size(v, 1);

v1 = B(1);
v2 = B(idivide(numel(B), int32(2)) + 1);

uv0 = conformal_nature_multi(v, f, B, H, 1:nv, v1, v2);

figure;
drawmesh(f, uv0, B);


