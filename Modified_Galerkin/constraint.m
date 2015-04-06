function [c,ceq] = constraint(x, L, lambda, epsilon)
% Contraint function for basis transformation

% Regenerate transformation matrix
X = x*(x'*x)^(-1/2);

% Produce new linear interation terms and lambda matrix
L_tilde = X'*L*X;
lambda_tilde = (X'*diag(lambda)*X);

% Set value of equality constraint
ceq = sum(sum(L_tilde.*lambda_tilde))- epsilon;
c = [];
end