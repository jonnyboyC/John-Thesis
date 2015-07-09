function [residual, X, C_til, L_til, Q_til, modal_amp_til, t] = ...
    optimal_rotation(epsilon, C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, evals)
%% Find Optimal Rotation 

timers = tic;
[X, flag] = constrained_POD(lambda, L, num_modes, basis_modes, epsilon, evals);
%inputs of constrained_POD are the POD temporal coefficients,modal_amp',
%the linear Galerkin matrix, l, the transformation dimensions N and
%n, and the transfer term parameter epsi
lambda_til = X'*diag(lambda(1:num_modes))*X;

% Modified reduced order model coefficients
L_til = X'*L*X;
C_til = X'*C;
Q_til=zeros(basis_modes,basis_modes,basis_modes);
Q = reshape(Q,num_modes,num_modes,num_modes);

for i = 1:basis_modes
    for j = 1:basis_modes
        for k = 1:basis_modes
            for p = 1:num_modes
                for r = 1:num_modes
                    for q = 1:num_modes
                        Q_til(i,j,k) = Q_til(i,j,k) + X(p,i)*Q(p,q,r)*X(q,j)*X(r,k);
                    end
                end
            end
        end
    end
end
Q_til=reshape(Q_til,basis_modes,basis_modes*basis_modes);

%% Numerically Find Error

% Galerkin 1D coefficients using term order
Gal_coeff=[C_til L_til Q_til]; 
Mod = true;

reduced_model_coeff = ode_coefficients2(basis_modes, Gal_coeff, Mod);

%Solution of the system of equation
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Calculate values at same time intervals 
 ao = til_init(modal_amp, init, basis_modes, X);
[t, modal_amp_til] = ode113(@(t, y) system_odes_mod2(t, y, reduced_model_coeff)...
    , tspan, ao, options); 

% move from non-dimension back to normal
t = t/t_scale;

% Calculate difference 
residual = error_til(lambda_til, modal_amp_til, t, tspan, flag);

% Report time
toc(timers);
end