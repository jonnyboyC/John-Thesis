function [integration, time] = time_integration(odesolver, system, integration, modal_TKE, i, ao, t_scale, tspan, options)
%#ok<*PFBNS>
tic;

vis = system{i}.vis;

m = flow_comps(system{i}.coef);
m_comps = flow_ncomps(system{i}.coef);

cnt = 1;

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

% Generate waiting bar
h = waitbar(0, 'Performing time integration on models in parrallel');

% Fetch results as they become available
for j = 1:cnt-1
    % fetch results if it take over half an hour discard results
    [~, t_job, modal_amp_job, type, subtype, time] = fetchNext(parfeval_futures);
    
    integration{i}.t.(type).(subtype) = t_job/t_scale;
    integration{i}.modal_amp.(type).(subtype) = modal_amp_job;
    
    % update wait bar
    waitbar(j/(cnt-1), h, sprintf('Galerkin system with %s %s finished in %6.1fs', ...
        strrep(type, '_', ' '), strrep(subtype, '_', ' '), time));
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
        system_odes_NL2(t, y, coef, eddy, vis, modal_TKE), ...
        tspan, ao, options);
    time = toc;
else
    tic;
    [t_job, modal_amp_job] = odesolver(@(t,y) ...
        system_odes2(t, y, coef), ...
        tspan, ao, options);
    time = toc;
end
end