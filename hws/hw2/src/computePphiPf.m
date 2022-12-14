function PphiPf = computePphiPf(F)
%% Compute partial derivative PphiPf of symmetric Dirichlet energy
%% Args:
%%      F[2, 2]: deformation graident
%% Returns:
%%      PphiPf[4, 1]: partial derivative

%% polar decomposition
[R, S] = polarDecomposition(F);

%% invariants
I1 = trace(S);
I2 = norm(S, "fro")^2;
I3 = det(S);

%% PPhiPI
PphiPI = [0; 0.5*(1+1/I3^2); -I2/I3^3];

%% PIPf
PIPf = zeros(4, 3);
PIPf(:, 1) = vec(R);
PIPf(:, 2) = 2*vec(F);

twist = [0, -1; 
         1,  0];
G = twist * F * twist';
PIPf(:, 3) = vec(G);

%% PphiPf
PphiPf = PIPf * PphiPI;

end
