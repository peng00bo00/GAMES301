function gamma = bestFitCurve(L, Ltarget, k)
%% Find the best fitted curve with given length and exterior angles
%% Args:
%%      L[nB, 1]: length of original curve segment
%%      Ltarget[nB, 1]: length of target curve segment
%%      k[nB, 1]: exterior angle of target curve segment
%% Returns:
%%      gamma[nB, 2]: fitted curve vertex

%% accumulate exterior angles and tangent vectors
phi = cumsum(k) - k(1);
T = [cos(phi) sin(phi)]';

%% boundary mass matrix
l = 0.5*(L + circshift(L, 1));     %% dual length
Ninv = diag(l);

%% solve optimal length to close the curve
TNTinv = matrixInv2x2(T * Ninv * T');
L = Ltarget - Ninv*T'*TNTinv*T*Ltarget;

%% accumulate scaled tangents
gamma = cumsum(L .* T', 1);
gamma = circshift(gamma, 1, 1);
gamma(1,:) = [0 0];

end
