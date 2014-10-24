function [rep] = error_til(Lam_til, a_til)

 rep = sum(diag(Lam_til))-mean(sum(a_til.^2,2));
end
 
 