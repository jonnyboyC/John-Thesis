function [X, L, epsilon, lambda] = constrained_POD(a,li_i,N,n,epsi)
L=li_i;
epsilon = epsi;
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
options.Algorithm = 'interior-point';
problem.options = options;

[x,~,~,OUTPUT,~] = fmincon(problem);
OUTPUT.message

X = x*(x'*x)^(-1/2);
end