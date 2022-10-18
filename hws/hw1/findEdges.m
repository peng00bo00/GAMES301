function E = findEdges(F)
%% find edges of a mesh
nF = size(F, 1);

E = zeros(3*nF, 2);

for i=1:nF
    f = sort(F(i,:));
    vi = f(1); vj = f(2); vk = f(3);
    
    E(3*i-2, :) = [vi, vj];
    E(3*i-1, :) = [vi, vk];
    E(3*i  , :) = [vj, vk];
end

%% remove repeated rows
E = unique(E,'rows');

end