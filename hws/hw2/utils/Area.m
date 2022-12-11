function A = Area(V)
%% Compute area of a given triangle
%% Args:
%%      V[3, 2]: each row is a vertex coordinates
%% Returns:
%%      A: area of the triangle
e1 = V(2, :) - V(1, :);
e2 = V(3, :) - V(1, :);

A = 0.5*norm(cross([e1 0], [e2 0]));

end