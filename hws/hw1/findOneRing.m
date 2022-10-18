function ring = findOneRing(heIdx, HE, NEXT, TWIN)
%% find one ring of vi
ring = [];

%% get one half-edge starting from vi
idx= heIdx;

%% loop over one-ring
while 1
    %% add vj to the list
    he = HE(idx, :);
    vj = he(2);
    ring = [ring, vj];
    
    %% move to next half-edge
    idx = NEXT(TWIN(idx));

    if idx == heIdx
        break
    end
end

end