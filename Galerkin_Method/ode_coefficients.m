function mod_coeff = ode_coefficients(num_modes, Gal_coeff, MOD)
% ODE_COEFFICENTS convert the system coefficeients into a single row for 
% use in time integration
%
%   mod_coeff = ODE_COEFFICIENTS(num_modes, Gal_coeff, MOD)


if MOD
    % Prefill array, copy linear terms, create offsets for quadratic terms
    mod_coeff = zeros(num_modes, num_modes+1 + (num_modes*num_modes));
    mod_coeff(:,1:num_modes+1) = Gal_coeff(:,1:num_modes+1);
    idx = num_modes + 2;
    offset = num_modes + 1;
else
    % Prefill array, copy linear terms, create offsets for quadratic terms
    mod_coeff = zeros(num_modes, num_modes + (num_modes*num_modes));
    mod_coeff(:,1:num_modes) = Gal_coeff(:,1:num_modes);
    idx = num_modes + 1;
    offset = num_modes;
end

% used by sub2ind
mat_size = [num_modes num_modes];

% Add terms to the mod_coeff array, sum terms mirrored across the diagonal
for i = 1:num_modes;
    for j = 1:num_modes;
        mod_coeff(:,idx) = Gal_coeff(:,offset + sub2ind(mat_size, j, i));
        idx = idx + 1;
    end
end
end

