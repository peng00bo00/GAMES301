function uv = BFFAuto(V, F, B)
%% Automatic parameterization with free boundary
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%%      B[1, nB]: boundary vertex index
%% Returns:
%%      uv[nV, 2]: uv coordinates

u = zeros(length(B), 1);
uv = BFFScale(V, F, B, u);

end