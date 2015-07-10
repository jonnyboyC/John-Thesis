function [rep] = error_til(Lam_til, a_til, t, tspan, flag)
% Total error in eliminating crtical transfer term
rep = sum(diag(Lam_til))-mean(sum(a_til.^2,2));

% pentalize for not fully integrating
if ~isequal(t, tspan)
    rep = (rep)*((length(tspan)/length(t))^2);
end

if flag == -2
    rep = rep*1e6;
end
end
 
 