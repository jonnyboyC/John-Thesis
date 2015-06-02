function mod_coeff = ode_coefficients(num_modes, Gal_coeff, MOD)
% Creating linear row of coefficients for eachs modes linear and quadratic
% terms
% ODE_COEFFICENTS(NUM_MODES, FCUH) created a row vector of coefficients for
% every mode in FCUH for NUM_MODES modes.

if MOD
    % Prefill array, copy linear terms, create offsets for quadratic terms
    mod_coeff = zeros(num_modes, num_modes+1 + (num_modes)*(num_modes+1)/2);
    mod_coeff(:,1:num_modes+1) = Gal_coeff(:,1:num_modes+1);
    idx = num_modes + 2;
    offset = num_modes + 1;
else
    % Prefill array, copy linear terms, create offsets for quadratic terms
    mod_coeff = zeros(num_modes, num_modes + (num_modes)*(num_modes+1)/2);
    mod_coeff(:,1:num_modes) = Gal_coeff(:,1:num_modes);
    idx = num_modes + 1;
    offset = num_modes;
end

% used by sub2ind
mat_size = [num_modes num_modes];

% Add terms to the mod_coeff array, sum terms mirrored across the diagonal
for i = 1:num_modes;
    for j = i:num_modes;
        if i == j
            mod_coeff(:,idx) = Gal_coeff(:,offset + sub2ind(mat_size, j, i));
        else
            mod_coeff(:,idx) = Gal_coeff(:,offset + sub2ind(mat_size, j, i)) ...
                             + Gal_coeff(:,offset + sub2ind(mat_size, i, j));
        end
        idx = idx + 1;
    end
end
end
