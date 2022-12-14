%% clear
clc; clear;

%% add path
addpath("./src");
addpath("./utils");

%% parameters
maxSteps= 100;
maxIter = 100;
lam     = 1e-5;
c       = 1e-5;
tau     = 0.9;

%% read mesh
meshName = 'ex';
path = fullfile('./mesh/', meshName);
[V, F] = readObj(path);
[B, H] = findBoundary(V, F);
nV = size(V, 1);

% %% drawing
figure;
h = drawmesh(F, V, B);

%% tutte's embedding
uv = tutte(V, F);

%% drawing
figure;
h = drawmesh(F, uv, B);

%% projected Newton solver
uv = projectedNewtonSolver(V, F, uv, maxSteps, maxIter, lam, c, tau);

% path = fullfile('./gif/', meshName);
% uv = projectedNewtonSolverSave(V, F, uv, maxSteps, maxIter, lam, c, gamma, path);

%% drawing
figure;
h = drawmesh(F, uv, B);

path = fullfile('./results/', [meshName '.obj']);
writeOBJ(path, [uv, zeros(nV, 1)], F);