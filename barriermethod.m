% Convert the objective into barrier form
t = orig_t;
while true
    % Create barrier function
    phi = -sum(log(-constraints));
    obj_phi = t * objective + phi;

    % Obtain its gradient and hessian
    obj_phi_grad = gradient(obj_phi, vars);
    obj_phi_hess = hessian(obj_phi, vars);
    
    while true
        % Substitute numbers
        grad_x = double(subs(obj_phi_grad, vars, x));
        hess_x = double(subs(obj_phi_hess, vars, x));
        
        % Calculate the Newton step (and the dual)
        obj_mat = [hess_x A'; A zeros(num_h, num_h)];
        val_mat = -[grad_x; zeros(num_h, 1)];
        
        res = linsolve(obj_mat, val_mat);
        newton_step = res(1:length(vars));
        newton_dual = res(length(vars)+1:length(res));
        
        % Calculate newton decrement (without the sqrt)
        % Note: instead of performing inverse on hessian the matrix
        % division is used for performance and accuracy reasons
        newton_dec = newton_step' * hess_x * newton_step;
        
        if newton_dec / 2 <= eps
            break
        end
        
        % Perform backtracking line search
        inner_t = 1;
        
        while true
            % While the search function value lies above the tangent line 
            % (of factor alpha), make t smaller
        	fn_val = subs(obj_phi, vars, (x' + inner_t * newton_step)');
            line_i = subs(obj_phi, vars, x) +...
                     alpha * inner_t * grad_x' * newton_step;

            if ~isreal(fn_val) || ~isreal(line_i) || fn_val > line_i
                inner_t = beta * inner_t;
            else
                break
            end
        end
        
        % Update x
        x = x + inner_t * newton_step';
    end
        
    % Check stopping criterion
    if m / t < eps
        break
    end
        
    % Increase t
    t = mu * t;
end