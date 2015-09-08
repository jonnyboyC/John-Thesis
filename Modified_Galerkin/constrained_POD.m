function [X, flag] = constrained_POD(lambda, L, n, N, epsilon, evals)
% CONSTRAINED_POD minimum rotation of the POD basis, attempt to produce a
% rotation that follows the energy constraints

% Initial guess for transformation matrix
x0 = zeros(N,n);
x0(1:n,1:n) = eye(n,n);

% Set up optimization problem
problem.objective = @(x) objective(x, lambda);
problem.nonlcon = @(x) constraint(x, L, lambda, epsilon);
problem.x0 = x0;
problem.solver = 'fmincon';

options = optimoptions('fmincon');
options.MaxFunEvals = evals;
options.Display = 'off';
problem.options = options;

% Perform minimization
[x,~,flag,OUTPUT,~] = fmincon(problem);
OUTPUT.message

% Calculate full transformation
X = x*(x'*x)^(-1/2);
end