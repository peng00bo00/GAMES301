function J = findJacobian(x1, x2)
%% Find Jacobian matrix from rest pose x1 to current pose x2
%% Args:
%%      x1[3, 2]: rest pose
%%      x2[3, 2]: current pose
%% Return:
%%      J[2, 2]: Jacobian matrix (deformation gradient)

Dm = (x1(2:end, :) - x1(1, :))';
Ds = (x2(2:end, :) - x2(1, :))';

%% J * Dm = Ds
J = Ds * matrixInv2x2(Dm);

end