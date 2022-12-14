function PfPx = computePfPx(x1)
%% Compute partial derivative of deformation gradient
%% Args:
%%      x1: rest pose
%% Returns:
%%      PfPx: partial derivative of deformation gradient

Dm = (x1(2:end, :) - x1(1, :))';
DmInv = matrixInv2x2(Dm);
a = DmInv(1, 1); c = DmInv(1, 2);
b = DmInv(2, 1); d = DmInv(2, 2);

PfPx = zeros(4, 6);

% pfpx0 = [-1,-1; 0, 0] * DmInv;
% pfpx1 = [ 0, 0;-1,-1] * DmInv;
% pfpx2 = [ 1, 0; 0, 0] * DmInv;
% pfpx3 = [ 0, 0; 1, 0] * DmInv;
% pfpx4 = [ 0, 1; 0, 0] * DmInv;
% pfpx5 = [ 0, 0; 0, 1] * DmInv;
% 
% PfPx(:, 1) = vec(pfpx0);
% PfPx(:, 2) = vec(pfpx1);
% PfPx(:, 3) = vec(pfpx2);
% PfPx(:, 4) = vec(pfpx3);
% PfPx(:, 5) = vec(pfpx4);
% PfPx(:, 6) = vec(pfpx5);

%% hard coded PfPx
PfPx(1, 1) = -(a+b); PfPx(1, 2) =      0; PfPx(1, 3) = a; PfPx(1, 4) = 0; PfPx(1, 5) = b; PfPx(1, 6) = 0;
PfPx(2, 1) =      0; PfPx(2, 2) = -(a+b); PfPx(2, 3) = 0; PfPx(2, 4) = a; PfPx(2, 5) = 0; PfPx(2, 6) = b;
PfPx(3, 1) = -(c+d); PfPx(3, 2) =      0; PfPx(3, 3) = c; PfPx(3, 4) = 0; PfPx(3, 5) = d; PfPx(3, 6) = 0;
PfPx(4, 1) =      0; PfPx(4, 2) = -(c+d); PfPx(4, 3) = 0; PfPx(4, 4) = c; PfPx(4, 5) = 0; PfPx(4, 6) = d;

end