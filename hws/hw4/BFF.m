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
A = A + 1e-5*speye(nV);

%% seperatate interior and boundary vertex
I = setdiff(1:nV, B);

AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

%% discrete curvatures
[K, kappa] = discreteCurvature(V, F, B);

%% target curvature at the boundary
k = zeros(nB, 1);

%% solve Yamabe equation (Neumann to Dirichlet)
phi = K; phi(B) = phi(B) - (kappa - k);
u = A \ phi;
u = u-u(1);     %% constant offset

%% solve boundary curve length
uB = u(B); 
uBr = circshift(uB, -1); Br = circshift(B, -1);     %% left shift the array to find right endpoint
L  = vecnorm(V(B,:)-V(Br,:), 2, 2);
Ltarget  = exp((uB+uBr)/2) .* L;

%% best fit curve
gamma = bestFitCurve(L, Ltarget, k);

%% extend curve
uv = zeros(nV, 2);

end