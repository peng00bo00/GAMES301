function [R, S, U, Sigma, V] = polarSVD(F)
%% Same as svd_rv, but returns polar decomposition as well
[U, Sigma, V] = svd_rv(F);

R = U * V';
S = V * Sigma * V';

end