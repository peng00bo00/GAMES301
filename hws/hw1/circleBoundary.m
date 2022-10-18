function VB = circleBoundary(B)
% Place the boundary points on a circle

nB = size(B, 2);
VB = zeros(nB, 2);

theta = (0:nB-1) * (2*pi/nB);

VB(:,1) = cos(theta);
VB(:,2) = sin(theta);

end