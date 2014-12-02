function [c,ceq] = constraint(x, L, lambda, epsilon)
% Contraint function for optimal rotation

 X = x*(x'*x)^(-1/2);
 L_tilde = X'*L*X;
 lambda_tilde = (X'*diag(lambda)*X);
 ceq = sum(sum(L_tilde.*lambda_tilde))- epsilon;
 c = [];
 end