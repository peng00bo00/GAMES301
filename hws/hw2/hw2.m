%% clear
clc; clear;

%% read mesh
[V, F] = readObj('cathead');
nV= size(V, 1);