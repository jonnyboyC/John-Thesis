function [c,ceq] = constraint(x)
 global L lambda epsilon;
 X = x*(x'*x)^(-1/2);
 L_tilde = X'*L*X;
 lambda_tilde = (X'*diag(lambda)*X);
 ceq = sum(sum(L_tilde.*lambda_tilde))- epsilon;
 c = [];
 end