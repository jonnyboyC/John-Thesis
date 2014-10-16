function dy = system_odes(t, y, model_coeff)
num_modes = (size(model_coeff,1));
num_terms = (size(model_coeff,2));
dy = zeros(num_modes, 1);

for i = 1:num_modes
    idx = 1;
    dy(i) = model_coeff(i,1);
    for j = 1:num_modes
        idx = idx + 1;
        dy(i) = dy(i) + model_coeff(i, idx)*y(j);
    end
    % index upper half of matrix
    for j = 1:num_modes
        for k = j:num_modes
            idx = idx + 1;
            dy(i) = dy(i) + model_coeff(i,idx)*y(j)*y(k);
        end
    end
end