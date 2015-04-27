function [rep] = error_til(Lam_til, a_til, t, tspan)
% Total error in eliminating crtical transfer term
rep = sum(diag(Lam_til))-mean(sum(a_til.^2,2));

% pentalize for not fully integrating
if ~isequal(t, tspan)
    rep = rep*10^(length(tspan)/length(t));
end
end
 
 