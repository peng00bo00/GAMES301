function b = computeGradient(V, F, uv, X1, As, PfPxs)
%% gradient of symmetric Dirichlet energy
%% Args:
%%      V[nV, 3]: vertex coordinates in 3D
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%% Returns:
%%      b[nV*2, 1]: gradient of current parameterization

nV = size(V, 1);
nF = size(F, 1);

b = zeros(nV, 2);

%% loop over all the triangles
for i=1:nF
    %% F[i, :] at local frame
    x1 = X1(:,:,i);
    q  = As(i);

    %% F[i, :] at uv plane
    x2= uv(F(i, :), :);
    
    %% deformation gradient from local frame to uv plane
    f = findJacobian(x1, x2);
    
    %% partial derivative of phi
    PphiPf = computePphiPf(f);
    PfPx   = PfPxs(:,:,i);

    PphiPx = q * PfPx' * PphiPf;

    v1 = F(i, 1); v2 = F(i, 2); v3 = F(i, 3);
    b(v1, 1) = b(v1, 1) + PphiPx(1); b(v1, 2) = b(v1, 2) + PphiPx(2);
    b(v2, 1) = b(v2, 1) + PphiPx(3); b(v2, 2) = b(v2, 2) + PphiPx(4);
    b(v3, 1) = b(v3, 1) + PphiPx(5); b(v3, 2) = b(v3, 2) + PphiPx(6);
end

b = reshape(b, [nV*2, 1]);

end