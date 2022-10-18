%% clear
clc; clear;

%% read mesh
[V, F] = readObj('cathead');
nV= size(V, 1);
B = findBoundary(V, F);

%% uniform weights
uv = tutte(V, F);

%% floater weights
uv = floater(V, F);

%% draw
h = drawmesh(F, uv);

% %% write to obj
% UV = zeros(nV, 3);
% UV(1:end, 1:2) = uv;
% writeOBJ("parametrization.obj", UV, F);