function reduced_mod_coeff = ode_coefficients(num_modes, desired_modes, fcuh)
% Loop to generate full model coefficients and selectively pick
% coefficients for use in a reduced model

mod_coeff = model_coefficients(num_modes, fcuh);

% Pull off reduced num_modes, from mod_coeff
%% TODO figure out why dm + 1
reduced_mod_coeff = zeros(desired_modes, (desired_modes^2)/2+ ceil(1.5*desired_modes)+1);
reduced_mod_coeff(1:desired_modes, 1:desired_modes +1) ...
    = mod_coeff(1:desired_modes, 1:desired_modes+1);


idx = desired_modes + 1;
for i = 1:desired_modes
    for j = i:desired_modes
        % index up half of matrix that is in column form
        idx = idx + 1;
        reduced_mod_coeff(:,idx) = mod_coeff(1:desired_modes, ...
                                i/2*(2*num_modes-i+3) + (j-i+1));
    end
end
end

function mod_coeff = model_coefficients(num_modes, fcuh)
% Prefill with zeros, pull first num_modes + 1 columns straight from fcuh
mod_coeff = zeros(size(fcuh,1), (num_modes^2)/2+ceil(1.5*num_modes)+1);
mod_coeff(:,1:num_modes+1) = fcuh(:,1:num_modes + 1);


idx = num_modes + 1;
offset = idx;
mat_size = [num_modes num_modes];

% Loop through remaining num_modes^2 columns. Sum if off diagonal terms so
% that mod_coeff = fcuh_ij + fcuh_ij
for i = 1:num_modes;
    for j = i:num_modes;
        idx = idx + 1;
        if i == j
            mod_coeff(:,idx) = fcuh(:,offset + sub2ind(mat_size, j, i));
        else
            mod_coeff(:,idx) = fcuh(:,offset + sub2ind(mat_size, j, i)) ...
                             + fcuh(:,offset + sub2ind(mat_size, i, j));
        end
    end
end
end

