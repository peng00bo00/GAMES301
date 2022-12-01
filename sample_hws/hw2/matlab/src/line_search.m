function [new_uv, step, success] = line_search(uv, d, min_step_to_singularity, wolfe_c1, wolfe_c2, e, initial)
    % input: e (handle of energy class)
    % output: uv_new

    max_iter = 50;
    min_step = 0.0;
    max_step = 0.99 * min_step_to_singularity;

    step = initial;

    old_e = e.compute_energy(uv);
    old_derivate = dot(d, e.compute_gradient(uv)); % < 0

    new_uv = uv + step * d;
    new_e = e.compute_energy(new_uv);

    for i=1:max_iter        
        if new_e > old_e + wolfe_c1 * step * old_derivate
            max_step = step;
            step = (min_step + max_step) / 2;
        else
            new_derivate = dot(d, e.compute_gradient(new_uv));

            if abs(new_derivate) > wolfe_c2 * abs(old_derivate)
                if new_derivate > 0
                    max_step = step;
                else
                    min_step = step;
                end
                step = (min_step + max_step) / 2;
            else
                success = true;
                return;
            end
        end

        new_uv = uv + step * d;
        new_e = e.compute_energy(new_uv);
    end

    success = false;
end

