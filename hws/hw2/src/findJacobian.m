function J = findJacobian(x1, x2)
%% find Jacobian matrix from rest pose x1 to new pose x2
%% each row of x1/x2 is the coordinate of the triangle vertex

Dm = (x1(2:end, :) - x1(1, :))';
Ds = (x2(2:end, :) - x2(1, :))';

%% J * Dm = Ds
J = Ds * matrixInv2x2(Dm);

end