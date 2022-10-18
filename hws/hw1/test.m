[x, t] = readObj('cathead');  % read triangle mesh from file

B = findBoundary(x, t);   % find the boundary loop of the mesh

uv = x(:, 2:3);   % trivial parameterization 
distortion = x(:, 1);  % arbitrary distortion for visulization


%% draw
B = findBoundary(x, t);   % find the boundary loop of the mesh
figure; subplot(121); h = drawmesh(t, x, B);
subplot(122); h2 = drawmesh(t, uv);

set(h2, 'CData', distortion, 'facecolor', 'interp');  % 
