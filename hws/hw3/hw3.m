%% clear
clc; clear;

%% read mesh
meshName = 'cow';
path = fullfile('./mesh/', meshName);
[V, F] = readObj(path);
[B, ~] = findBoundary(V, F);

nV = size(V, 1);

%% LSCM
uv1 = LSCM(V, F);

path = fullfile('./results/LSCM', [meshName '.obj']);
writeOBJ(path, [uv1, zeros(nV, 1)], F);

%% Mean Value Coordinates
uv2 = MVCTutte(V, F);

path = fullfile('./results/MVC', [meshName '.obj']);
writeOBJ(path, [uv2, zeros(nV, 1)], F);

figure;
drawmesh(F, uv1, B);
 
figure;
drawmesh(F, uv2, B);