function gamma = bestFitCurve(L, Ltarget, k)
%% Find the best fitted curve with given length and exterior angles
%% Args:
%%      L[nB, 1]: length of original curve segment
%%      Ltarget[nB, 1]: length of target curve segment
%%      k[nB, 1]: exterior angle of target curve segment
%% Returns:
%%      gamma[nB, 2]: fitted curve vertex

nB = length(L);

N = diag(1./L);
Ninv = diag(L);

gamma = zeros(nB, 2);

end