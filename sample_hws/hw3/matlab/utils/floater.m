function uv = floater(v, f, B)

% input: v (vertices), f (triangles, mx3), B (boundary indices)
% output: parametrization for each vertex

nv = size(v, 1);
in = setdiff(1:nv, B);

%% Get VV-F matrix and part of half-edge data structure
VV2F = VertexVertexToFace(v, f);
VV2N = VertexVertexToNext(VV2F, f);

VV2Ft = VV2F';

%% for every interior vertex, calculate floater weights

Mat = sparse(nv, nv);
Vec = zeros(nv, 2);

uv_bound = get_coords_bound(v, B);
bounds = zeros(nv, 2);
bounds(B, :) = uv_bound;

is_bound = zeros(nv, 1);
is_bound(B) = 1;

weiMat = sparse(nv, nv); 

for i = 1:numel(in)
    
    x = in(i); % don't use i below
    Y = find(VV2Ft(:, x), 1);  % neighbour of x, but not ordered
    lst = [Y(1)];   % neighbour of x, ordered
    y = lst(1);
    
    while true
        y = VV2N(y, x);
        if y == lst(1)
            break
        end
        
        lst(end + 1) = y;
    end
    
    wei = floater_weights(v(x, :), v(lst, :));
    weiMat(x, lst) = wei;
    
    for id = 1:numel(lst)
       y = lst(id);
       if (is_bound(y))
           Vec(x, :) = Vec(x, :) + wei(id) * bounds(y, :);
       else
           Mat(y, x) = -wei(id);
       end
    end
end

%% construct and solve linear system

Mat = Mat';
Mat = Mat(in, in);
Mat = Mat + speye(numel(in));

Vec = Vec(in, :);

uv_in = [Mat \ Vec(:, 1), Mat \ Vec(:, 2)];

uv = zeros(nv, 2);
uv(in, :) = uv_in;
uv(B, :) = uv_bound;

% check

% tmp_uv = weiMat * uv;
% 
% for x = 1:nv
%     if (~is_bound(x))
%         tmp_uv(x, :) - uv(x, :)
%     end
% end

end

%% Vertex-Vertex to Face Matrix
function VV2F = VertexVertexToFace(v, f)
    nv = size(v, 1);
    [nf, polyDim] = size(f);
    VV2F = sparse(f, f(:, [2:end 1]), repmat(1:nf, 1, polyDim), nv, nv);
end

%% Vertex-Vertex to Next Vertex Matrix

function k = findNext(j, fid, f)
    polyDim = size(f, 2);

    idx = find(arrayfun(@(x) x == j, f(fid, :)), 1);
    idx = mod(idx, polyDim) + 1;
    k = f(fid, idx);
end

function VV2N = VertexVertexToNext(VV2F, f)
    nv = size(VV2F, 1);
    VV2N = sparse(nv, nv);
   
    for j = 1:nv
        I = find(VV2F(:, j));
        for i_id = 1:numel(I)
            i = I(i_id);
            VV2N(i, j) = findNext(j, VV2F(i, j), f);
        end
    end
end

%% Get FLoater Weights

function A = floater_weights(o, x)

% input: x (center of subgraph), xi (neigbour of x)
% output: w (local weights, s.t. x = \sum_{i) w(i) * xi(i))

n = size(x, 1);

if n <= 2
    warning("number of neighbours <= 2");
end

for i = 1:n
    x(i, :) = x(i, :) - o;
end

length = arrayfun(@(i) norm(x(i, :)), 1:n);
angle = zeros(1, n);

for i = 1:n
    j = mod(i, n) + 1;
    angle(j) = angleVV(x(i, :), x(j, :)); % angle between x(i), x(j)
end

angle = 2 * pi * angle / sum(angle);

for i = 2:numel(angle)
    angle(i) = angle(i - 1) + angle(i);
end

ps = [reshape(length.*cos(angle), n, 1), reshape(length.*sin(angle), n, 1)];

M = zeros(n, n);

for i = 1:n
    for j = 1:n
        M(i, j) = cross2(ps(i, :), ps(j, :));
    end
end

A = zeros(n, 1);

for i = 1:n
    for j = 1:n
        k = mod(j, n) + 1;
        if (j ~= i && k ~= i)
            if (M(i, j) >= 0 && M(i, k) < 0) % [i, j, k] include O
                sk = M(i, j); % k
                si = M(j, k); % i
                sj = M(k, i); % j
                S = sk + si + sj;
                A([i, j, k]) = A([i, j, k]) + [si; sj; sk] / S;
                break
            end
        end
    end
end

A = A / n;

% check:
% all = [0, 0];
% for i = 1:n
%     all = all + A(i) * ps(i, :);
% end

end

function ang = angleVV(x, y)
    ang = atan2(norm(cross(x, y)), dot(x, y));
end

function v = cross2(a, b)
    v = a(1) * b(2) - a(2) * b(1);
end

%% parametrization of boundary vertices, circle for an example
function uv = get_coords_bound(v, B)
    % input: v (set of vertices), B (list of boundary point)
    % output: uv (size(B) x 2, coords of boundary)
    
    n = numel(B);
    uv = reshape([cos(2*pi*(1:n)'/n), sin(2*pi*(1:n)'/n)], n, 2);
    
end

function A = uniform_weights(o, x)
    n = size(x, 1);
    A = ones(n, 1) / n;
end
