function uv = BFFAuto(V, F)
%% Automatic parameterization with free boundary
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      uv[nV, 2]: uv coordinates

%% find boundary
[B, ~] = findBoundary(V, F);

u = zeros(length(B), 1);
uv = BFFScale(V, F, u);

end