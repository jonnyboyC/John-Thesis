function X = constrained_POD(lambda, L, N, n, epsilon, evals)
% Function to determine the optimal rotation of the system

x0 = eye(n,n);
if N > n
    x0(N,:) = 0;
end

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
%options.Algorithm = 'active-set';


[x,~,~,OUTPUT,~] = fmincon(problem);
OUTPUT.message

X = x*(x'*x)^(-1/2);
end