function [R, S] = polarDecomposition(F)
%% polar decomposition with rotation-variant SVD
%% copy from Dynamic Deformables: Implementation and Production Practicalities
%% http://www.tkim.graphics/DYNAMIC_DEFORMABLES/

%% rotation variant SVD
[U, Sigma, V] = svd_rv(F);

R = U * V';
S = V * Sigma * V';

end

function [U, Sigma, V] = svd_rv(F)
%% rotation-variant SVD
%% copy from Dynamic Deformables: Implementation and Production Practicalities
%% http://www.tkim.graphics/DYNAMIC_DEFORMABLES/

%% standard SVD
[U, Sigma, V] = svd(F);

%% reflection matrix
%% modified for 2D deformation
L = eye(2,2);
L(2,2) = det(U * V');

%% see where to pull the reflection out of
detU = det(U);
detV = det(V);

if (detU < 0 && detV > 0)
    U = U * L;
elseif (detU > 0 && detV < 0)
    V = V * L;
end

%% push the reflection to the diagonal
Sigma = Sigma * L;

end