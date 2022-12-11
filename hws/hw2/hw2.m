%% clear
clc; clear;

%% add path
addpath("./src");
addpath("./utils");

%% parameters
maxSteps= 3000;
maxIter = 100;
lam     = 1e-5;
c       = 1e-4;
gamma   = 0.9;

%% read mesh
[V, F] = readObj('./mesh/camelhead.obj');
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
uv = projectedNewtonSolver(V, F, uv, maxSteps, maxIter, lam, c, gamma);
% uv = projectedNewtonSolverSave(V, F, uv, 1000, maxIter, lam, c, gamma, './gif/hand');

%% drawing
figure;
h = drawmesh(F, uv, B);

% writeOBJ("./results/cow.obj", [uv, zeros(nV, 1)], F);