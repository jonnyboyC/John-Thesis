function [t, modal_amp] = time_integration(reduced_model_coeff, eddy, modal_TKE, i, t, modal_amp, options)
%#ok<*PFBNS>

ao = modal_amp_flux(init,1:num_modes)+modal_amp_mean(init, 1:num_modes);
t_temp = cell(size(t,1),1);
modal_amp_temp = cell(size(modal_amp,1),1);
num_models = size(eddy,1);

% Integrate Galerkin Systems of requested methods using linear eddy visocity
parfor j = 1:num_models
    if ~isempty(reduced_model_coeff{j,1});
        fprintf('Performing ode113 on Galerkin system with %s with linear eddy visocity\n', reduced_model_coeff{j,2}); 

        tic;
        [t_temp{j}, modal_amp_temp{j}] = ode113(@(t,y) ...
            system_odes(t, y, reduced_model_coeff{j,1}, eddy{j,1,i}, Re0, modal_TKE, true), ...
            tspan, ao, options);  
        toc2 = toc;

        fprintf('Completed in %f6.4 seconds\n\n', toc2);
        modal_amp_temp{j} = modal_amp_temp{j}(:,2:end);
    end
end


% Integrate Galerkin Systems of requested methods using nonlinear eddy
% visocity
parfor j = 3:num_models
    if ~isempty(reduced_model_coeff{j,1});
        fprintf('Performing ode113 on Galerkin system with %s with nonlinear eddy visocity\n', reduced_model_coeff{j,2}); 

        tic;
        [t_temp{j+num_models}, modal_amp_temp{j+num_models}] = ode113(@(t,y) ...
            system_odes_NL(t, y, reduced_model_coeff{j,1}, eddy{j,1,i}, Re0, modal_TKE, true), ...
            tspan, ao, options);  
        toc2 = toc;

        fprintf('Completed in %f6.4 seconds\n\n', toc2);
        modal_amp_temp{j+num_models} = modal_amp_temp{j+num_models}(:,2:end);
    end
end

t(:,:,i) = [t_temp, reduced_model_coeff(:,2)];
modal_amp(:,:,i)  = [modal_amp_temp, reduced_model_coeff(:,2)];