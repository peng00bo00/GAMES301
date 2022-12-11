function H = projectHessian(V, F, uv, X1, As, PfPxs)
%% Projected Hessian of symmetric Dirichlet energy
%% Args:
%%      V[nV, 3]: vertex coordinates in 3D
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%% Returns:
%%      H[nV*2, nV*2]: Hessian matrix of current parameterization

nV = size(V, 1);
nF = size(F, 1);

%% sparse matrix indices and values
I = zeros(nF*36, 1); J = zeros(nF*36, 1);
Value = zeros(nF*36, 1);

%% update with each triangle
for i=1:nF
    %% F[i, :] at local frame
    x1 = X1(:,:,i);
    q  = As(i);

    %% F[i, :] at uv plane
    x2= uv(F(i, :), :);
    
    %% deformation gradient from local frame to uv plane
    f = findJacobian(x1, x2);

    %% SVD
    [eigs, lams] = SDEigenSystem(f);
    eigs = max(eigs, 0);

    %% local Hessian
%     Hq = zeros(4);
%     for ii = 1:4
%         Hq = Hq + eigs(ii) * lams(:, ii) * lams(:, ii)';
%     end

    Hq = lams * diag(eigs) * lams';

    %% PfPx
    PfPx   = PfPxs(:,:,i);
    Hq = q * PfPx' * Hq * PfPx;

    %% update indices and values
    v1x = F(i, 1); v2x = F(i, 2); v3x = F(i, 3);
    v1y = v1x+nV;  v2y = v2x+nV;  v3y = v3x+nV;

    indices = [v1x, v1y, v2x, v2y, v3x, v3y];
    idxI = zeros(6) + indices';
    idxJ = zeros(6) + indices;
    
    idx = (i-1)*36+1;   % starting index
    I(idx:idx+35)     = reshape(idxI, 1, []);
    J(idx:idx+35)     = reshape(idxJ, 1, []);
    Value(idx:idx+35) = reshape(Hq,   1, []);

end

H = sparse(I, J, Value, nV*2, nV*2);

end