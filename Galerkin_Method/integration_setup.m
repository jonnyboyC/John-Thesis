function [reduced_model_coeff] = integration_setup(eddy, vis, l, q, idx, total_models, linear_models, num_modes)

reduced_model_coeff = cell(total_models,2);
Gal_coeff   = cell(total_models,2);
Mod = false;

% Final setup for time integration
for j = 1:total_models
    if ~isempty(eddy{j,1,idx})

        % Setup up system coefficients
        if j <= linear_models
            total_vis = repmat((eddy{j,1,idx}+vis), 1, size(l{j,1,idx},1));
            Gal_coeff{j,1,idx} = [l{j,1,idx}.*total_vis, q{j,1,idx}];
            Gal_coeff{j,2,idx} = eddy{j,2,idx};
        else
            Gal_coeff{j,1,idx} = [l{j,1,idx}, q{j,1,idx}];
            Gal_coeff{j,2,idx} = eddy{j,2,idx};
        end

        % Galerkin 1D coefficients using 
        reduced_model_coeff{j,1} = ode_coefficients(num_modes, Gal_coeff{j,1,idx}, Mod);
        reduced_model_coeff{j,2} = Gal_coeff{j,2,idx};
    end
end
    
end