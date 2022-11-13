clc;

%% initialize 2 triangles
% x1 = [[0,-1];[1,0];[0,1]];
% x2 = [[1, 2];[0,2];[0,-1]];
x1 = rand(3,2);
x2 = rand(3,2);

%% deformation gradient
F = findJacobian(x1, x2);
f = vec(F);


i=3;j=2;
delta = 1e-7;

xx2 = x2;
xx2(i,  j) = xx2(i, j)+delta;
FF = findJacobian(x1, xx2);
ff1 = vec(FF);

xx2 = x2;
xx2(i,  j) = xx2(i, j)-delta;
FF = findJacobian(x1, xx2);
ff2 = vec(FF);


(ff1 - ff2) ./ (2*delta)

ComputePfPx(x1, x2)