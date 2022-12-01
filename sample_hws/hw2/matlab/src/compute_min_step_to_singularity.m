function min_step_to_singularity = compute_min_step_to_singularity(D, X, dx)

Dx = full(D * X);
Ddx = full(D * dx);

nf = size(Dx, 1) / 4;

J_11 = Dx(1:nf);         % J_11
J_21 = Dx(nf+1:2*nf);    % J_21
J_12 = Dx(2*nf+1:3*nf);  % J_12
J_22 = Dx(3*nf+1:4*nf);  % J_22

dJ_11 = Ddx(1:nf);           % dJ_11
dJ_21 = Ddx(nf+1:2*nf);      % dJ_21
dJ_12 = Ddx(2*nf+1:3*nf);    % dJ_12
dJ_22 = Ddx(3*nf+1:4*nf);    % dJ_22

a = dJ_11 .* dJ_22 - dJ_21 .* dJ_12;
b = dJ_11.*J_22 - dJ_12.*J_21 + dJ_22.*J_11 - dJ_21.*J_12;
c = J_11.*J_22 - J_12.*J_21;

[min_pos_rt, ~, ~] = quadratic_equation(a, b, c);
min_step_to_singularity = max(min(min_pos_rt), 0);

if (imag(min_step_to_singularity) ~= 0)
    fprintf('\n\twarning: complex???\n\n');
    min_step_to_singularity
end

end

% solve ax^2 + bx + c = 0, x1 < x2
function [min_pos_rt, x1, x2] = quadratic_equation(a, b, c)
    x1 = zeros(size(a));
    x2 = zeros(size(a));
    min_pos_rt = zeros(size(a));
    min_pos_rt(:) = Inf;

    sqdelta = sqrt(b.^2 - 4 * (a .* c));

%     eps = 1e-13;
%     ill = find(abs(a) < eps | abs(b) < eps | abs(c) < eps);
%     if(numel(ill) > 0)
%         fprintf('\twarning: %d ill triangle\n', numel(ill));
% %         for i=1:numel(ill)
% %             fprintf('\twarning: ill triangle: a: %.15f, b: %.15f, c: %.15f\n', a(ill(i)), b(ill(i)), c(ill(i)));
% %         end               
%     end

    pos = a > 1e-15;
    neg = a < -1e-15;
    zero = ~pos & ~neg;

    x1(pos) = (-b(pos) - sqdelta(pos)) ./ (2 * a(pos));
    x2(pos) = (-b(pos) + sqdelta(pos)) ./ (2 * a(pos));

    x1(neg) = (-b(neg) + sqdelta(neg)) ./ (2 * a(neg));
    x2(neg) = (-b(neg) - sqdelta(neg)) ./ (2 * a(neg));

    x1(zero & (b ~= 0)) = -c(zero & (b ~= 0)) ./ b(zero & (b ~= 0));
    x1(zero & (b == 0)) = Inf;
    x2(zero) = Inf;

    
    x2_pos = imag(x2) == 0 & real(x2) > 0;
    x1_pos = imag(x1) == 0 & real(x1) > 0;
    min_pos_rt(x2_pos) = x2(x2_pos);
    min_pos_rt(x1_pos) = x1(x1_pos);
end

