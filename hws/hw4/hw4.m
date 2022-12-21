clc; clear;

%% read mesh
meshName = 'cathead';
path = fullfile('./mesh/', meshName);

[V, F] = readObj(path);
[B, ~] = findBoundary(V, F);

nV = size(V, 1);
nF = size(F, 1);
nB = length(B);

I = setdiff(1:nV, B);