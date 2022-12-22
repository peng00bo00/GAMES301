function uv = BFFScale(V, F, B, u)
%% Boundary First Flattening with given scale factors
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%%      B[1, nB]: boundary vertex index
%%      u[nB, 1]: target scale factor at the boundary
%% Returns:
%%      uv[nV, 2]: uv coordinates

nV = size(V, 1);

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

%% solve Yamabe equation (Dirichlet to Neumann)
phi =-K;
a = AII \ (phi(I)-AIB*u);
k = kappa - (phi(B) - AIB'*a - ABB*u);
k = k / sum(k) * 2 * pi;

%% solve boundary curve length
uB = u;
uBr = circshift(uB, -1); Br = circshift(B, -1);     %% left shift the array to find right endpoint
L  = vecnorm(V(B,:)-V(Br,:), 2, 2);
Ltarget  = exp((uB+uBr)/2) .* L;

%% best fit curve
gamma = bestFitCurve(L, Ltarget, k);

%% extend curve
uv = zeros(nV, 2);
uv(B, 1) = gamma(:, 1);
uv(I, 1) = AII \ (-AIB*gamma(:, 1));

%% Hilbert transform
aB = uv(B, 1); h = zeros(nV, 1);
h(B, :) =-0.5*(circshift(aB, -1) - circshift(aB, 1));
uv(:, 2) = A \ h;

uv = uv - mean(uv, 1);

end