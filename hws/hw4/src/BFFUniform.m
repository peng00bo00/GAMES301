function uv = BFFUniform(V, F)
%% Uniformization to a disk with BFF
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      uv[nV, 2]: uv coordinates

%% find boundary
[B, ~] = findBoundary(V, F);

nV = size(V, 1);
Br = circshift(B, -1);

for i=1:10
    %% dual edge length
    L = vecnorm(V(B,:)-V(Br,:), 2, 2);
    l = 0.5*(L + circshift(L, 1));

    %% set exterior angle proportional to the most recent dual lengths
    k = l / sum(l) * 2 * pi;
    
    uv = BFFAngle(V, F, k);
    V = [uv zeros(nV, 1)];
end

uv = uv(:, 1:2);

end