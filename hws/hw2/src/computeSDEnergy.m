function E = computeSDEnergy(x1, x2)
%% Symmetric Dirichlet energy from x1 to x2
%% Args:
%%      x1: rest pose
%%      x2: current pose
%% Return:
%%      E: symmetric Dirichlet energy

F = findJacobian(x1, x2);
Finv = matrixInv2x2(F);

E = norm(F, "fro")^2 + norm(Finv, "fro")^2;
E = 0.5*E;

end