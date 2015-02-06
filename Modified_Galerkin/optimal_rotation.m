function [rep, X, C_til, L_til, Q_til] = ...
    optimal_rotation(epsilon, ci, li, q, OG_nm, RD_nm, lambda2, eig_func, t, init, evals)
% TODO update comments   

tic;
X = constrained_POD(eig_func',li,OG_nm,RD_nm,epsilon,evals);
%inputs of constrained_POD are the POD temporal coefficients,eig_func',
%the linear Galerkin matrix, li, the transformation dimensions N and
%n, and the transfer term parameter epsi
Lam_til = X'*diag(lambda2(1:OG_nm))*X;

% Modified reduced order model coefficients
L_til = X'*li*X;
C_til = X'*ci;
Q_til=zeros(RD_nm,RD_nm,RD_nm);
Q = reshape(q,OG_nm,OG_nm,OG_nm);
for i = 1:RD_nm
    for j = 1:RD_nm
        for k = 1:RD_nm
            for p = 1:OG_nm
                for r = 1:OG_nm
                    for q = 1:OG_nm
                        Q_til(i,j,k) = Q_til(i,j,k) + X(p,i)*Q(p,q,r)*X(q,j)*X(r,k);
                    end
                end
            end
        end
    end
end
% solve new system of equation
Q_til=reshape(Q_til,RD_nm,RD_nm*RD_nm);
fc_til=[C_til L_til Q_til]; % constant, linear and quadratic terms coefficients

reduced_model_coeff= ode_coefficients(RD_nm,RD_nm,fc_til);

%Solution of the system of equation
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
%(Control implicit in the modes)

% Look more into tspan;
tspan = t;
[~,Y_til] = ode113(@(t, y) system_odes(t, y, -reduced_model_coeff)...
    , tspan ,eig_func(init,1:RD_nm), options); 

rep = error_til(Lam_til,Y_til);
toc;
end