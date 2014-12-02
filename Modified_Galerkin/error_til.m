function [rep] = error_til(Lam_til, a_til)
% Total error in eliminating crtical transfer term
 rep = sum(diag(Lam_til))-mean(sum(a_til.^2,2));
end
 
 