function [Fsub, I, B, H] = make_holes(V, F, x, num)

nv = size(V, 1);
nf = size(F, 1);
B = findBoundary(V, F);
I = setdiff(1:nv, B);

[Vout, ~, Fout, ~] = local_patch(V, F, x, num);

Fin = ~Fout;
Fsub = F(Fin, :);

[B, H] = findBoundary(V, Fsub);
I = setdiff(setdiff(1:nv, Vout), B);
for i=1:numel(H)
    I = setdiff(I, H{i});
end

end

