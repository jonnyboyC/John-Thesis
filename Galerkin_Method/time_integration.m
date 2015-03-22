function [t, modal_amp] = time_integration(reduced_model_coeff, eddy, Re0, modal_TKE, i, t, modal_amp, ao, tspan, total_models, linear_models, options)
%#ok<*PFBNS>

t_temp = cell(size(t,1),1);
modal_amp_temp = cell(size(modal_amp,1),1);

% Integrate Galerkin System usign async parallel execution
parfeval_futures = parallel.FevalOnAllFuture;
for j = 1:total_models
    parfeval_futures(j) = parfeval(@integration, 3, reduced_model_coeff{j,1}, ...
        eddy{j,1,i}, Re0, modal_TKE, ao, tspan, j, linear_models, options);
end

% Generate waiting bar
h = waitbar(0, 'Performing time integration on models in parrallel');

% Fetch results as they become available
for j = 1:total_models
    % fetch results if it take over half an hour discard results
    [job_idx, t_job, modal_amp_job, time] = fetchNext(parfeval_futures, 3000);
    
    t_temp(job_idx) = {t_job};
    modal_amp_temp(job_idx) = {modal_amp_job};
    
    % update wait bar
    waitbar(j/total_models, h, sprintf('Galerkin system with %s finished in %6.1fs', reduced_model_coeff{j,2}, time));
end

% In case any don't finish fill with zeros
for j = 1:total_models
    if isempty(t_temp{j})
       t_temp{j} = 0;
       modal_amp_temp = zeros(size(ao));
    end
end

close(h);
t(:,:,i) = [t_temp, reduced_model_coeff(:,2)];
modal_amp(:,:,i)  = [modal_amp_temp, reduced_model_coeff(:,2)];

end

function [t_job, modal_amp_job, time] = integration(reduced_model_coeff, eddy, Re0, modal_TKE, ao, tspan, j, linear_models, options)
if ~isempty(reduced_model_coeff); 
    if j <= linear_models
        % Time integration with linear eddy visocity
        tic;
        [t_job, modal_amp_job] = ode113(@(t,y) ...
            system_odes(t, y, reduced_model_coeff), ...
            tspan, ao, options);  
        time = toc;
    else
        % Time inetgration with nonlinear eddy visocity
        tic;
        [t_job, modal_amp_job] = ode113(@(t,y) ...
            system_odes_NL(t, y, reduced_model_coeff, eddy, Re0, modal_TKE), ...
            tspan, ao, options);  
        time = toc;  
    end
    modal_amp_job = modal_amp_job(:,2:end);
end
end