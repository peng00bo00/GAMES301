function [Xb] = setBoundary(n, shape)
    if strcmp(shape, 'circle') % generate circle boundary
        theta = ((1:n) .* (2 * pi / n)).';
        Xb = [cos(theta), sin(theta)];

    elseif strcmp(shape, 'rec') % generate rectangle boundary
        m = double(idivide(int32(n), 4));   % \floor(n, 4)
        r = mod(n, 4);
        num = ones(1, 4) .* m;
        num(1:r) = num(1:r) + 1;            % number of points on every edge of square 

        cons = @(v, k) (ones(1,k).*v).';
        lins = @(a, b, k) (linspace(a, b - (b - a) / k, k).'); % uniform in [a, b)

        A = [lins(-1, 1, num(1)), cons(-1, num(1))]; % bottom
        B = [cons(1, num(2)), lins(-1, 1, num(2))];  % right
        C = [lins(1, -1, num(3)), cons(1, num(3))];  % top
        D = [cons(-1, num(4)), lins(1, -1, num(4))]; % left

        Xb = cat(1, cat(1, cat(1, A, B), C), D);

    else
        warning('invalid parameter of shape, should be ''circle'' or ''rec''');
    end
end