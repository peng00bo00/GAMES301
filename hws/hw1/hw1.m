%% clear
clc; clear;

%% read mesh
[V, F] = readObj('cathead');
nV= size(V, 1);
B = findBoundary(V, F);

%% uniform weights
uv = tutte(V, F);

%% draw
h1 = drawmesh(F, uv);

% %% write to obj
% VTutte = zeros(nV, 3);
% VTutte(1:end, 1:2) = uv;
% writeOBJ("tutte.obj", VTutte, F)

%% floater weights
uv = floater(V, F);

%% draw
h2 = drawmesh(F, uv);

% %% write to obj
% VFloater = zeros(nV, 3);
% VFloater(1:end, 1:2) = uv;
% writeOBJ("floater.obj", VFloater, F)