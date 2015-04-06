function goal = objective(x, lambda)
% Generate a new matrix from the pseudo inverse to garentee that X'X = I
X = x*(x'*x)^(-1/2);

% Take the sum of lambdas original ie original energy
sum_lambda = sum(lambda);

% Take the new energy sum from the transformed lambdas
sum_lambda_tilde = sum(diag(X'*diag(lambda)*X));

% Minimize energy differnce
goal = (sum_lambda - sum_lambda_tilde)/sum_lambda;
end