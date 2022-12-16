%% clear
clc; clear;

%% read mesh
meshName = 'cathead';
path = fullfile('./mesh/', meshName);
[V, F] = readObj(path);
[B, ~] = findBoundary(V, F);

%% LSCM
uv = LSCM(V, F);

figure;
drawmesh(F, uv, B);