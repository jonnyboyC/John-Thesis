function [X, L, lambda] = constrained_POD(a,L,N,n,epsilon)
lambda = zeros(N, 1);
for i=1:N;
    lambda(i) = mean(a(i,:).*a(i,:));
end

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
% options.Algorithm = 'active-set';
% options.MaxFunEvals = 30000;
problem.options = options;

[x,~,~,OUTPUT,~] = fmincon(problem);
OUTPUT.message

X = x*(x'*x)^(-1/2);
end