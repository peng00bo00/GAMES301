function [W, R] = local_step(J)

% input:  J (Jacobian for each face [d, d, f])
% output: W (weight for each face   [d, d, f])
%         R (local step result      [d, d, f])

[U, Sj, V] = pagesvd(J, "vector");
Sw = sqrt((1 + Sj.^(-1)).*(1 + Sj.^(-2)));

% fprintf("max_sig: %.4f\n", max(Sw, [], 3));

W = pagemtimes(U, (Sw .* pagetranspose(U)));
R = pagemtimes(U, "none", V, "transpose");

end