filename = 'cathead.obj';
file_big = 'camelhead.obj';
file_2d = '2D.obj';
[v, f, ~, ~] = readObj(file_big, false);
[B, H] = findBoundary(v, f);

tic;
uv = floater(v, f, B);
% uv = tutte(v, f, B);
toc;

figure;
drawmesh(f, uv, B);


