function [R, S] = polarDecomposition(F)
%% polar decomposition with rotation-variant SVD
%% copy from Dynamic Deformables: Implementation and Production Practicalities
%% http://www.tkim.graphics/DYNAMIC_DEFORMABLES/

%% rotation variant SVD
[U, Sigma, V] = svd_rv(F);

R = U * V';
S = V * Sigma * V';

end
