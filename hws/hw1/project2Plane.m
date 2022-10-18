function [P, radius, Theta] = project2Plane(V, vi, ring)
%% project vertices to the plane with geodesic polar map
nP = length(ring);
P  = zeros(nP, 2);

Theta = zeros(nP+1, 1);
radius= zeros(nP, 1);

%% place the first point
radius(1) = norm(V(ring(1), :) - V(vi,:));

%% loop over other points
for vj=2:nP
%     e1 = V(ring(vj-1), :) - V(vi,:); e1 = e1 / norm(e1);
%     e2 = V(ring(vj), :)   - V(vi,:); e2 = e2 / norm(e2);
% 
%     theta(vj) = theta(vj-1) + acos(dot(e1, e2));
    Theta(vj) = Theta(vj-1) + findAngle(V(vi, :), V(ring(vj), :), V(ring(vj-1), :));
    radius(vj) = norm(V(ring(vj), :) - V(vi,:));
end

% e1 = V(ring(1),   :) - V(vi,:); e1 = e1 / norm(e1);
% e2 = V(ring(end), :) - V(vi,:); e2 = e2 / norm(e2);
% theta(end) = theta(end-1) + acos(dot(e1, e2));
Theta(end) = Theta(end-1) + findAngle(V(vi, :), V(ring(1), :), V(ring(end), :));

%% normalize theta
Theta = Theta / Theta(end) * (2*pi);
Theta = Theta(1:end-1, :);

P(:, 1) = radius .* cos(Theta);
P(:, 2) = radius .* sin(Theta);

end

function theta = findAngle(vi, vj, vk)
%% a helper function to find angle vj-vi-vk
e1 = vj - vi; e1 = e1 / norm(e1);
e2 = vk - vi; e2 = e2 / norm(e2);

theta = acos(dot(e1, e2));

end