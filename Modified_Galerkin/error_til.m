function residual = error_til(lambda_til, a_til, t, tspan, flag)
% ERROR_TIL calculate the energy from numerical intergration, penalize
% trasnfromed models that either don't integrate the whole time length or
% are from invalid rotation
%
%   residual = ERROR_TIL(lambda_til, a_til, t, tspan, flag)

% Total error in eliminating crtical transfer term
residual = sum(diag(lambda_til))-mean(sum(a_til.^2,2));

% pentalize for not fully integrating
if ~isequal(t, tspan)
    residual = (residual)*((length(tspan)/length(t))^2);
end

% penalize for coming from an invalid rotation
if flag == -2
    residual = residual*1e6;
end
end
 
 