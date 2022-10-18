function [bj, bk, bl] = findBaryCententricCoordinate(pj, pk, pl)
%% find barycentric coordinate of the origin on the triangle
ejk = [pk - pj, 0]; ejl = [pl - pj, 0];
eij = [pj, 0]; eik = [pk, 0]; eil = [pl, 0];

%% actually this is the doubled area
A = norm(cross(ejk, ejl));

Aj = norm(cross(eik, eil));
Ak = norm(cross(eij, eil));
Al = norm(cross(eij, eik));

bj = Aj/A; bk = Ak/A; bl = Al/A;

end