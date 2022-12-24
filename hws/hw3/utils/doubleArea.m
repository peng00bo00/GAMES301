function As = doubleArea(V, F)
%% Compute doubled area of triangles
%% Args:
%%      V[nV, 3]: vertices in 3D
%%      F[nF, 3]: face connectivity
%% Returns:
%%      As[nF, 1]: doubled area of triangles

nF = size(F, 1);

Es = reshape(V(F(:, 2:3), :)-V(F(:, [1 1]), :), [nF 2 3]);
As = vecnorm(cross(Es(:,1,:), Es(:, 2, :), 3), 2, 3);

end