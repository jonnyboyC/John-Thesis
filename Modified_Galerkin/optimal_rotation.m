function [residual, X, C_til, L_til, Q_til, modal_amp_til, t] = ...
    optimal_rotation(epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, evals)
%% Find Optimal Rotation 

timers = tic;
[X, flag] = constrained_POD(lambda, L, OG_nm, RD_nm, epsilon, evals);
%inputs of constrained_POD are the POD temporal coefficients,modal_amp',
%the linear Galerkin matrix, l, the transformation dimensions N and
%n, and the transfer term parameter epsi
lambda_til = X'*diag(lambda(1:OG_nm))*X;

% Modified reduced order model coefficients
L_til = X'*L*X;
C_til = X'*C;
Q_til=zeros(RD_nm,RD_nm,RD_nm);
Q = reshape(Q,OG_nm,OG_nm,OG_nm);

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
Q_til=reshape(Q_til,RD_nm,RD_nm*RD_nm);

%% Numerically Find Error

% Galerkin 1D coefficients using term order
Gal_coeff=[C_til L_til Q_til]; 
Mod = true;

reduced_model_coeff = ode_coefficients2(RD_nm, Gal_coeff, Mod);

%Solution of the system of equation
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Calculate values at same time intervals 
 ao = til_init(modal_amp, init, RD_nm, X);
[t, modal_amp_til] = ode113(@(t, y) system_odes_mod2(t, y, reduced_model_coeff)...
    , tspan, ao, options); 

% Calculate difference 
residual = error_til(lambda_til, modal_amp_til, t, tspan, flag);

% Report time
toc(timers);
end