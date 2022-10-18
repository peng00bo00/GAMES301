function [HE, NEXT, TWIN, VHE] = findHE(V, F)
%% find half-edges of a mesh
nV = size(V, 1);
nF = size(F, 1);

nHE  = 3*nF;
HE   = zeros(nHE, 2);
NEXT = zeros(nHE, 1);

%% vertex pointer to half-edges
VHE  =-ones(nV, 1);

%% add half-edges for each face
for i=1:nF
    vi = F(i, 1); vj = F(i, 2); vk = F(i, 3);
    
    %% half-edges
    HE(3*i-2, 1) = vi; HE(3*i-2, 2) = vj;
    HE(3*i-1, 1) = vj; HE(3*i-1, 2) = vk;
    HE(3*i  , 1) = vk; HE(3*i  , 2) = vi;

    %% next HE
    NEXT(3*i-2) = 3*i-1;
    NEXT(3*i-1) = 3*i;
    NEXT(3*i  ) = 3*i-2;

    if VHE(vi) < 0
        VHE(vi) = 3*i-2;
    end

    if VHE(vj) < 0
        VHE(vj) = 3*i-1;
    end

    if VHE(vk) < 0
        VHE(vk) = 3*i;
    end
end

%% twin of HE
twin = zeros(nHE, 2);
twin(:, 1) = HE(:, 2);
twin(:, 2) = HE(:, 1);

[~, TWIN] = ismember(HE, twin, "rows");

end