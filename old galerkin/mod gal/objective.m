function [f] = objective(x)
 global lambda;
 X = x*(x'*x)^(-1/2);
 sum_lambda = sum(lambda);
 lambda_tilde = diag(X'*diag(lambda)*X);
 sum_lambda_tilde = sum(lambda_tilde);

 f = (sum_lambda - sum_lambda_tilde)/sum_lambda;
 end