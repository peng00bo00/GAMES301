function uv = BFFSquare(V, F, B)
%% Uniformization to a square with BFF
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%%      B[1, nB]: boundary vertex index
%%      k[nB, 1]: target curvature at the boundary
%% Returns:
%%      uv[nV, 2]: uv coordinates

nV = size(V, 1);

%% fix boundary exterior angle
nB= length(B);
k = zeros(nB,1);
k(floor((1:4)*nB/4)) = 2*pi/4;

%% build Laplacian matrix
A = cotLaplacian(V, F);
A = A + 1e-8*speye(nV);

%% seperatate interior and boundary vertex
I = setdiff(1:nV, B);

AII = A(I,I);
AIB = A(I,B);
ABB = A(B,B);

%% discrete curvatures
[K, kappa] = discreteCurvature(V, F, B);

%% solve Yamabe equation (Neumann to Dirichlet)
phi = -K;
h = kappa - k;      %% Neumann data
phi(B) = phi(B) - h;

u = A \ phi;
u = u-mean(u);      %% constant offset

%% solve boundary curve length
uB = u(B); 
uBr= circshift(uB, -1); Br = circshift(B, -1);     %% left shift the array to find right endpoint
L  = vecnorm(V(B,:)-V(Br,:), 2, 2);
Ltarget  = exp((uB+uBr)/2) .* L;

%% best fit curve
gamma = bestFitCurve(L, Ltarget, k);

%% extend curve with harmonics
uv = zeros(nV, 2);
uv(B, :) = gamma;
uv(I, :) = AII \ (-AIB*gamma);

uv = uv - mean(uv, 1);

end