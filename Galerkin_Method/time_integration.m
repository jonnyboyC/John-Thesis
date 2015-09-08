function [integration, time] = time_integration(odesolver, system, integration, modal_TKE, ...
                                            i, ao, t_scale, tspan, int_time, options)
% TIME_INTEGRATION time intergration of all POD-Galerkin models,
% intergration performed in parrallel 
%
%   [integration, time] = time_intergration(odesolver, system,
%       intergration, modal_TKE, i, ao, t_scale, tspan, int_time, options)
                                        
                                        
%#ok<*PFBNS>
tic;
pool = gcp('nocreate');

vis = system{i}.vis;        

m = flow_comps(system{i}.coef);
m_comps = flow_ncomps(system{i}.coef);

cnt = 1;
models = flow_fields(system{i}.eddy);

% Generate waiting bar
h = waitbar(0, 'Performing time integration on models in parrallel');

if ~isempty(pool)

    % Integrate Galerkin System usign async parallel execution
    parfeval_futures = parallel.FevalOnAllFuture;
    for j = 1:m_comps
        s = flow_comps(system{i}.coef.(m{j}));
        s_comps = flow_ncomps(system{i}.coef.(m{j}));

        for k = 1:s_comps
            parfeval_futures(cnt) = parfeval(@integrator, 5, odesolver, system{i}.coef.(m{j}).(s{k}), ...
                system{i}.eddy.(m{j}).(s{k}), vis, m{j}, s{k}, modal_TKE, ao, tspan, options);
            cnt = cnt + 1;
        end
    end

    % Fetch results as they become available
    for j = 1:models
        % fetch results if it take over an hour discard results
        [~, t_job, modal_amp_job, type, subtype, time] = fetchNext(parfeval_futures, int_time);

        if ~isempty(type) && ~isempty(subtype) && ~isempty(time)
            integration{i}.t.(type).(subtype) = t_job/t_scale;
            integration{i}.modal_amp.(type).(subtype) = modal_amp_job;

            % update wait bar
            waitbar(j/(models), h, sprintf('Galerkin system with %s %s finished in %6.1fs', ...
                strrep(type, '_', ' '), strrep(subtype, '_', ' '), time));
        end
    end
else
    for j = 1:m_comps
        s = flow_comps(system{i}.coef.(m{j}));
        s_comps = flow_ncomps(system{i}.coef.(m{j}));
        for k = 1:s_comps
            % Perform integration
            [t_job, modal_amp_job, type, subtype, time] = integrator(odesolver, ...
                system{i}.coef.(m{j}).(s{k}),system{i}.eddy.(m{j}).(s{k}), vis, ...
                m{j}, s{k}, modal_TKE, ao, tspan, options);
            
            if ~isempty(type) && ~isempty(subtype) && ~isempty(time)
                integration{i}.t.(type).(subtype) = t_job/t_scale;
                integration{i}.modal_amp.(type).(subtype) = modal_amp_job;

                % update wait bar
                waitbar(cnt/(models), h, sprintf('Galerkin system with %s %s finished in %6.1fs', ...
                    strrep(type, '_', ' '), strrep(subtype, '_', ' '), time));
            end
            cnt = cnt + 1;
        end
    end
end

cnt = 1;

% Fill in model/submodels with a structure telling integration took too long
for j = 1:m_comps
    s = flow_comps(system{i}.coef.(m{j}));
    s_comps = flow_ncomps(system{i}.coef.(m{j}));
    
    % Fill in blank model
    if ~isfield(integration{1}.t, m{j})
        integration{1}.t.(m{j}) = struct;
        integration{1}.modal_amp.(m{j}) = struct;
    end
    
    for k = 1:s_comps
        % Fill in blank submodels
        if ~isfield(integration{1}.t.(m{j}), s{k})
            integration{1}.t.(m{j}).(s{k}) = [];
            integration{1}.modal_amp.(m{j}).(s{k}) = [];
            cancel(parfeval_futures(cnt))
        end
        % Report any errors and add the error structure
        if ~isempty(pool) && ~isempty(parfeval_futures(cnt).Error)
            integration{i}.t.(m{j}).(s{k}).error = true;
            integration{i}.modal_amp.(m{j}).(s{k}).error = true;
            disp(parfeval_futures(cnt).Error.message);
        end
        cnt = cnt + 1;
    end
end
close(h);

time = toc;
end

function [t_job, modal_amp_job, type, subtype, time] = integrator(odesolver, coef, ...
                            eddy, vis, type, subtype, modal_TKE, ao, tspan, options)
                        
% If if a nonlinear model passed as an input
if strncmpi(subtype, 'NL', 2) 
    tic;
    [t_job, modal_amp_job] = odesolver(@(t,y) ...
        system_odes_NL(t, y, coef, eddy, vis, modal_TKE), ...
        tspan, ao, options);
    time = toc;
else
    tic;
    [t_job, modal_amp_job] = odesolver(@(t,y) ...
        system_odes(t, y, coef), ...
        tspan, ao, options);
    time = toc;
end
end