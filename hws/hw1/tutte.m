function uv = tutte(V, F)
%% Tutte's embedding with uniform weights
nV = size(V, 1);
nF = size(F, 1);

%% find boundary
B = findBoundary(V, F);
nB= size(B, 2);
VB= circleBoundary(B);

%% find edges
E = findEdges(F);
nE= size(E, 1);

%% uniform Laplacian L=D-A
I = []; J = []; value = [];
D = zeros(1, nV);

for e = 1:nE
    vi = E(e, 1); vj = E(e, 2);

    if ~ismember(vi, B)
        I = [I, vi];
        J = [J, vj];
        value = [value, -1];

        D(vi) = D(vi)+1;
    end

    if ~ismember(vj, B)
        I = [I, vj];
        J = [J, vi];
        value = [value, -1];

        D(vj) = D(vj)+1;
    end
end

%% add diagnoal element
I = [I, 1:nV];
J = [J, 1:nV];
value = [value, D];

%% add boundary constraints
I = [I, B];
J = [J, B];
value = [value, ones(1,nB)];

%% LHS
L = sparse(I, J, value, nV, nV);

%% RHS
b = zeros(nV, 2);
b(B, :) = VB;

%% solve linear system
uv= L \ b;

end