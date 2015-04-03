function mod_coeff = ode_coefficients(num_modes, Gal_coeff)
% Creating linear row of coefficients for eachs modes linear and quadratic
% terms
% ODE_COEFFICENTS(NUM_MODES, FCUH) created a row vector of coefficients for
% every mode in FCUH for NUM_MODES modes.

% Prefill array
mod_coeff = zeros(num_modes, num_modes + (num_modes)*(num_modes+1)/2);

% copy linear terms
mod_coeff(:,1:num_modes) = Gal_coeff(:,1:num_modes);

% Create offset for quadratic terms
idx = num_modes + 1;
offset = num_modes;

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

