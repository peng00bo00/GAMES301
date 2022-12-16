function uv = project2Plane(V)
%% Project a triangle to plane
%% Args:
%%      V[3, 3]: vertex coordinates in 3D, each row is a coordinate vector
%% Returns:
%%      uv[3, 2]: vertex coordinates in 2D

%% center in 3D
c = mean(V, 1);

%% normal
e1 = V(2, :) - V(1, :);
e2 = V(3, :) - V(1, :);
n = cross(e1, e2);
n = n / norm(n);

%% rotation
z = [0, 0, 1];
r = vrrotvec(n, z);

%% project to plane
uv = vrrotvec2mat(r) * (V - c)';
uv = uv(1:2, :)';

end