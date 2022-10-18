function uv = floater(V, F)
%% Tutte's embedding with Floater weights
nV = size(V, 1);
nF = size(F, 1);

%% find boundary
B = findBoundary(V, F);
nB= size(B, 2);
VB= circleBoundary(B);

%% find half-edges
[HE, NEXT, TWIN, VHE] = findHE(V, F);

%% Floater weighted Laplacian L
I = []; J = []; value = [];
D = zeros(1, nV);

for vi=1:nV
    if ~ismember(vi, B)
        %% find one half-edge starting from vi
        heIdx = VHE(vi);

        %% find one-ring
        ring = findOneRing(heIdx, HE, NEXT, TWIN);
        nRing= length(ring);

        %% project to plane
        [P, ~, Theta] = project2Plane(V, vi, ring);

        %% add weights
        for pjIdx=1:nRing
            %% find a valid triangle on the plane
            [pjIdx, pkIdx, plIdx] = find2DTriangle(pjIdx, Theta);
            vj = ring(pjIdx); vk = ring(pkIdx); vl = ring(plIdx);

            %% find barycentric coordinates
            pj = P(pjIdx, :); pk = P(pkIdx, :); pl = P(plIdx, :);
            [bj, bk, bl] = findBaryCententricCoordinate(pj, pk, pl);

            %% add to sparse matrix
            I     = [    I, vi, vi, vi];
            J     = [    J, vj, vk, vl];
            value = [value,-bj,-bk,-bl];

            D(vi) = D(vi)+1;
        end
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