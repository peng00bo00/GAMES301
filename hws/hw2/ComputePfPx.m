function PfPx = ComputePfPx(x1, x2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Dm = (x1(2:end, :) - x1(1, :))';
DmInv = inv(Dm);

PfPx = zeros(4, 6);

pfpx0 = [-1,-1; 0, 0] * DmInv;
pfpx1 = [ 0, 0;-1,-1] * DmInv;
pfpx2 = [ 1, 0; 0, 0] * DmInv;
pfpx3 = [ 0, 0; 1, 0] * DmInv;
pfpx4 = [ 0, 1; 0, 0] * DmInv;
pfpx5 = [ 0, 0; 0, 1] * DmInv;

PfPx(:, 1) = vec(pfpx0);
PfPx(:, 2) = vec(pfpx1);
PfPx(:, 3) = vec(pfpx2);
PfPx(:, 4) = vec(pfpx3);
PfPx(:, 5) = vec(pfpx4);
PfPx(:, 6) = vec(pfpx5);

end