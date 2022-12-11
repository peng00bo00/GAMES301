function [a] = vec(A)
%% vec operator on 2nd-order tensor
%% copy from Dynamic Deformables: Implementation and Production Practicalities
%% http://www.tkim.graphics/DYNAMIC_DEFORMABLES/

[rows, cols] = size(A);
a = reshape(A, rows * cols , 1);
end
