function [X, flag] = constrained_POD(lambda, L, N, n, epsilon, evals)
% Optimization problem to create an optimal rotation that minimizing the
% energy difference between the original system and the new system.

% Initial guess for transformation matrix
x0 = zeros(N,n);
x0(1:n,1:n) = eye(n,n);

% set up problem
problem.objective = @(x) objective(x, lambda);
problem.nonlcon = @(x) constraint(x, L, lambda, epsilon);
problem.x0 = x0;
problem.solver = 'fmincon';

options = optimoptions('fmincon');
options.MaxFunEvals = evals;
problem.options = options;

% Other options to consider changing
%options.UseParallel = true;
% options.Algorithm = 'sqp';

[x,~,flag,OUTPUT,~] = fmincon(problem);
OUTPUT.message

X = x*(x'*x)^(-1/2);
end