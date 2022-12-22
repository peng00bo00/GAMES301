clc; clear;

%% add path
addpath('./src');
addpath('./util');

%% read mesh
meshName = 'cow';
path = fullfile('./mesh/', meshName);

[V, F] = readObj(path);
[B, ~] = findBoundary(V, F);

nV = size(V, 1);
nF = size(F, 1);
nB = length(B);

% figure;
% drawmesh(F, V, B);

%% BFF
% uv = BFFUniform(V, F, B);
% uv = BFFAuto(V, F, B);
uv = BFFSquare(V, F, B);

figure;
drawmesh(F, uv, B);