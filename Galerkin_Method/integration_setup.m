function [reduced_model_coeff] = integration_setup(niu, l, q, i, num_models, num_modes)

reduced_model_coeff = cell(num_models,2);
Gal_coeff   = cell(num_models,2);

% Final setup for time integration
for j = 1:num_models
    if ~isempty(niu{j,1,i})

        % Setup up system coefficients
        Gal_coeff{j,1,i} = [l{j,1,i}, q{j,1,i}];
        Gal_coeff{j,2,i} = niu{j,2,i};

        % rearrage into 1D form
        reduced_model_coeff{j,1} = ode_coefficients(num_modes, Gal_coeff{j,1,i});
        reduced_model_coeff{j,2} = Gal_coeff{j,2,i};
    end
end
    
end