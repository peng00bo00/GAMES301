classdef SymDirichlet
    % class to compute energy and gradient for Sym Dirichlet Energy
    
    properties
        nv
        nf
        area          % area for each face
        grad_operator % matrix from uv to jacobian [J_11, J_12; J_21, J_22]
    end
    
    methods
        function obj = SymDirichlet(V, F)
            obj.nv = size(V, 1);
            obj.nf = size(F, 1);
            [obj.grad_operator, ~, obj.area] = deform_grad_operator(V, F);
        end

        function [jacobian] = compute_jacobian(self, uv)
            % output: jacobian [2 2 nf], [J_11, J_12; J_21, J_22] for each
            % block

            jacobian_u = full(self.grad_operator * sparse(uv(1 : self.nv)));           % J_11 J_12 u w.r.t xy
            jacobian_v = full(self.grad_operator * sparse(uv(self.nv+1 : 2*self.nv))); % J_21 J_22 v w.r.t xy
            jacobian = reshape([jacobian_u.'; jacobian_v.'], [2 2 self.nf]);
        end

        function [e] = compute_energy(self, uv)
            % require size(uv) == [nf * 2, 1]
            % output: e (energy) = 1/2 * \sum_f A_f (norm2(J)^2 + norm2(invJ)^2)

            jacobian = compute_jacobian(self, uv);
            inv_jacobian = matrix_inverse_22(jacobian);
            e = 0.5 * sum((jacobian.^2 + inv_jacobian.^2) .* reshape(self.area, [1 1 self.nf]), 'all');               
        end

        function [grad] = compute_gradient(self, uv)            
            jacobian = compute_jacobian(self, uv);
            inv_jacobian = matrix_inverse_22(jacobian);

            % transpose inv jacobian
            tinv_jacobian = inv_jacobian;
            tinv_jacobian(1, 2, :) = inv_jacobian(2, 1, :);
            tinv_jacobian(2, 1, :) = inv_jacobian(1, 2, :);

            % grad w.r.t jacobian
            grad_j = (jacobian - ...
                pagemtimes(pagemtimes(tinv_jacobian, inv_jacobian), tinv_jacobian));
            grad_j = grad_j .* reshape(self.area, [1 1 self.nf]);
            grad_j = reshape(grad_j, [2, 2 * self.nf]);

            % grad w.r.t. uv
            grad = reshape(full(sparse(grad_j) * self.grad_operator).', [self.nv * 2, 1]);
        end

        function [hessian] = compute_hessian(self, uv)
            % svd
            jacobian = compute_jacobian(self, uv);
            [U, S, V] = pagesvd(jacobian);
            
            % invarient
            I2 = S(1, 1, :).^2 + S(2, 2, :).^2;
            I3_inv = (S(1, 1, :) .* S(2, 2, :)).^(-1);
            
            % eigen value
            lam = zeros(1, 1, self.nf, 4);
            lam(:, :, :, 1) = 1 + 3 * S(1, 1, :).^(-4);
            lam(:, :, :, 2) = 1 + 3 * S(2, 2, :).^(-4);
            lam(:, :, :, 3) = 1 + I3_inv.^2 + I2 .* I3_inv.^3;
            lam(:, :, :, 4) = 1 + I3_inv.^2 - I2 .* I3_inv.^3;
            
            % eigen matrix
            ev = zeros(2, 2, self.nf, 4);
            ev(1, 1, :, 1) = 1;
            ev(2, 2, :, 2) = 1;
            ev(1, 2, :, 3) = 1;
            ev(2, 1, :, 3) = 1;
            ev(1, 2, :, 4) = -1;
            ev(2, 1, :, 4) = 1;

            for i=1:4
                ev(:, :, :, i) = pagemtimes(pagemtimes(U, ev(:, :, :, i)), 'none', V, 'transpose'); 
            end
            ev(:, :, :, 3:4) = ev(:, :, :, 3:4) .* sqrt(0.5);
            ev = reshape(ev, [4, 1, self.nf, 4]);
            
            % hessian of energy w.r.t F
            H_f = sum(pagemtimes(max(lam, 0) .* ev, 'none', ev, 'transpose'), 4); 
%             H_f = sum(pagemtimes(lam .* ev, 'none', ev, 'transpose'), 4); 
            weighted_H_f = reshape(self.area, [1 1 self.nf]) .* H_f;
            
            tripJ = reshape(repmat(1:4*self.nf, [4 1]), [4 4 self.nf]);
            tripI = pagetranspose(tripJ);

            weighted_H_diag = sparse( ...
                reshape(tripI, [self.nf*16, 1]), ...
                reshape(tripJ, [self.nf*16, 1]), ...
                reshape(weighted_H_f, [self.nf*16, 1]), self.nf*4, self.nf*4);
            
            % hessian of energy w.r.t X
            p = zeros(1, self.nf * 4);
            p(1 : 2 : self.nf*4) = 1 : self.nf*2;
            p(2 : 2 : self.nf*4) = self.nf*2+1 : self.nf*4;
            
            D = [self.grad_operator, sparse(self.nf * 2, self.nv);
                sparse(self.nf * 2, self.nv), self.grad_operator];

            hessian = D(p, :).' * weighted_H_diag * D(p, :);
        end
    end
end

function [invA] = matrix_inverse_22(A)
    % input: A   (shape = [2 2 n], block of 2x2 matrix)
    % output invA (shape = [2 2 n], each block contains inverse matrix)

    n = size(A, 3);
    detA = A(1, 1, :) .* A(2, 2, :) - A(1, 2, :) .* A(2, 1, :); % [1 1 n]
    
    invA = zeros(2, 2, n);
    invA(1, 1, :) = A(2, 2, :);
    invA(2, 2, :) = A(1, 1, :);
    invA(1, 2, :) = -A(1, 2, :);
    invA(2, 1, :) = -A(2, 1, :);

    invA = invA ./ detA;
end

