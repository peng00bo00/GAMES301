function uv = projectedNewtonSolverSave(V, F, uv, maxSteps, maxIter, lam, c, gamma, path)
%% One step projected Newton solver
%% Args:
%%      V[nV, 3]: vertex coordinates in 3D
%%      F[nF, 3]: face connectivity
%%      uv[nV, 2]: vertex coordinates in 2D
%%      maxSteps: maximum steps in projected Newton
%%      maxIter: maximum iteration in line search
%%      lam: positive-definite parameter
%%      c: search control parameters
%%      gamma: shrinkage factor in line search
%% Returns:
%%      uv_new[nV, 2]: updated uv coordinates

[B, H] = findBoundary(V, F);

nV = size(V, 1);
nF = size(F, 1);

%% reusable parameters
As    = zeros(nF, 1);          %% triangle area
X1    = zeros(3, 2, nF);    %% triangle at local frame
PfPxs = zeros(4, 6, nF);    %% PfPx at different triangles

for i=1:nF  
    x1 = project2Plane(V(F(i, :), :));

    As(i)        = Area(x1);
    X1(:,:,i)    = x1;
    PfPxs(:,:,i) = computePfPx(x1);
end

figure;
h = drawmesh(F, uv, B);
fname = fullfile(path, 'frame0');
saveas(gcf, fname, 'jpg')
close(gcf);

%% projected Newton step
for i=1:maxSteps
    fprintf('\nstep %d\n', i);
    
%     if stop
%         break
%     end
    
    tic;
    [uv, stop] = projectedNewtonStep(V, F, uv, X1, As, PfPxs, maxIter, lam, c, gamma);
    toc;
    
    if mod(i, 10) == 0
    figure;
    h = drawmesh(F, uv, B);
    fname = fullfile(path, ['frame' num2str(i)]);
    saveas(gcf, fname, 'jpg')
    close(gcf);
    end

end

end