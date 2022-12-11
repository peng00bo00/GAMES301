function E = computeTotalSDEnergy(F, uv, X1, As)
%% Compute total symmetric Dirichlet energy of the mesh with current parameterization
nF = size(F, 1);
E  = 0;

for i=1:nF
    %% F[i, :] at local frame
    x1 = X1(:,:,i);
    q  = As(i);
    
    %% F[i, :] at uv plane
    x2= uv(F(i, :), :);

    Eq = computeSDEnergy(x1, x2);
    E = E + q*Eq;
end

end