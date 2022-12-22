function Ainv = matrixInv2x2(A)
%% 2x2 matrix inverse
Ainv = zeros(2);
detA = A(1,1)*A(2,2) - A(1,2)*A(2,1);

Ainv(1,1) = A(2,2); Ainv(1,2) =-A(1,2);
Ainv(2,1) =-A(2,1); Ainv(2,2) = A(1,1);

Ainv = Ainv ./ detA;
end