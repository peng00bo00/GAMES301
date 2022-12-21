function uv = BFF(V, F, B)
%% Boundary First Flattening
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%%      B[1, nB]: boundary vertex index
%% Returns:
%%      uv[nV, 2]: uv coordinates

nV = size(V, 1);
nF = size(F, 1);
nB = length(B);

%% build Laplacian matrix
A = cotLaplacian(V, F);

%% seperatate interior and boundary vertex
I = setdiff(1:nV, B);

AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

%% solve Yamabe equation

%% best fit curve

%% extend curve

uv = zeros(nV, 2);

end