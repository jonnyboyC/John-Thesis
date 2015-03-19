function [t, modal_amp] = time_integration(reduced_model_coeff, eddy, Re0, modal_TKE, i, t, modal_amp, ao, tspan, total_models, linear_models, options)
%#ok<*PFBNS>

t_temp = cell(size(t,1),1);
modal_amp_temp = cell(size(modal_amp,1),1);

% Integrate Galerkin Systems of requested methods using linear eddy visocity
for j = 1:total_models
    if ~isempty(reduced_model_coeff{j,1});

        if j <= linear_models
            fprintf(['Performing ode113 on Galerkin system with %s',...
                'with linear eddy visocity\n'], reduced_model_coeff{j,2}); 
            tic;
            [t_temp{j}, modal_amp_temp{j}] = ode113(@(t,y) ...
                system_odes(t, y, reduced_model_coeff{j,1}), ...
                tspan, ao, options);  
            toc2 = toc;
        else
            fprintf(['Performing ode113 on Galerkin system with %s',...
                'with nonlinear eddy visocity\n'], reduced_model_coeff{j,2}); 
            tic;
            [t_temp{j}, modal_amp_temp{j}] = ode113(@(t,y) ...
                system_odes_NL(t, y, reduced_model_coeff{j,1}, eddy{j,1,i}, Re0, modal_TKE), ...
                tspan, ao, options);  
            toc2 = toc;  
        end

        fprintf('Completed in %f6.4 seconds\n\n', toc2);
        modal_amp_temp{j} = modal_amp_temp{j}(:,2:end);
    end
end

t(:,:,i) = [t_temp, reduced_model_coeff(:,2)];
modal_amp(:,:,i)  = [modal_amp_temp, reduced_model_coeff(:,2)];