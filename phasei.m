% Clear existing variables and output window
clear
clc

% Objective and constraints
syms x1 x2 x3

p2_objective = x1^4 + 10*x1^2 + 5*x2^2 + x3^2 + x1 + x2 + x3;
p2_constraints = [x1-5 3-x2 -x3];
orig_vars = symvar(p2_objective);

p2_A = [1 1 0; 0 1 5];
b = [10 20];

% Parameters
alpha = 0.2;
beta = 0.9;
mu = 10;
orig_t = 5;
eps = 0.0001;
m = length(p2_constraints);

% Starting points (must satisfy eq constraints; doesn't have to satisfy
% ineq constraints
sp = [[0 10 2]; [2 8 2]; [8 2 18/5]];
num_h = height(p2_A);
nu = zeros(1, num_h);

% Determine if the inequality consraint is satisfied

for idx = 1:length(sp)
    orig_x = sp(idx, :);
    s_val = double(max(subs(p2_constraints, orig_vars, orig_x)));
    disp("==========================")
    disp(['Current x: [' num2str(double(orig_x)) ']'])

    % If s < 0, the point is already strictly feasible, and no further 
    % action are needed. Otherwise, construct an an optimization problem 
    % that finds x  that makes s less than zero (i.e. phase I)
    if s_val >= 0
        disp("Current x value is not feasible, starting phase I...")
        syms s;
        p1_objective = s;
        p1_constraints = p2_constraints - s;

        objective = p1_objective;
        constraints = p1_constraints;
        A = [p2_A zeros((num_h), 1)];

        vars = [orig_vars s];
        x = [orig_x s_val+0.1];
        barriermethod
        
        disp(['The updated x value is: [' num2str(double(x)) ']'])
    end
        
    % Run the barrier method (phase II)
    disp("Current x value is already feasible, starting phase II...")
    objective = p2_objective;
    constraints = p2_constraints;
    A = p2_A;

    vars = orig_vars;
    if s_val >= 0 
        x = x(1:length(x)-1);
    else
        x = orig_x;
    end

    barriermethod
    disp(['Final x: [' num2str(double(x)) ']'])
    disp(['f0: ' num2str(double(subs(objective, vars, x))) ' f0_phi: ' num2str(double(subs(obj_phi, vars, x)))])
end