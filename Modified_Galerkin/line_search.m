function [epsilon_low, epsilon_high, transfer, found_flip] = line_search(problem)
% LINE_SEARCH brute force sweep of an area around the initial condition
% suggested by balajawecz.
%
% [epsilon_low, epsilon_high, transfer, flip] = LINE_SEARCH(problem) edit
% line_search for input details

% Unpack structure
C           = problem.C;
L           = problem.L;
Q           = problem.Q;
num_modes   = problem.num_modes;
basis_modes = problem.basis_modes;
lambda      = problem.lambda;
modal_amp   = problem.modal_amp;
t_scale     = problem.t_scale;
tspan       = problem.tspan;
init        = problem.init;
line_range  = problem.line_range;


% intial guess for epsilon
epsilon_0 = sum(L(1:num_modes, 1:num_modes)*lambda(1:num_modes)); 
epsilon = zeros(line_range+1,1);
transfer = zeros(line_range+1,1);
epsilon(1) = epsilon_0;
found_flip = false;

% Find initial value for line search and plot
transfer(1) = optimal_rotation(epsilon_0, C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, 64000);

% Brute force Line Search, rough pass parrallel loop. Idea is to quickly
% sweep a large area to look for sign changes to use for a finer pass in
% fzero

% Find out if a pool is currently created
pool = gcp('nocreate');
if ~isempty(pool)
    
    % Create jobs for the points in the vicinity of our initial guess
    parfeval_futures = parallel.FevalOnAllFuture;
    for i = 1:line_range
        if mod(i,2) 
            flip = 1;
        else
            flip = -1;
        end
        amp = floor(i+1/2)/128;
        epsilon(i+1) = epsilon_0*(1-amp*flip);
        parfeval_futures(i) = parfeval(@optimal_rotation, 1, epsilon(i+1), ...
            C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, 6000);
    end


    for i = 1:line_range
        
        % Retreive results as they become available
        [job_idx, transfer_job] = fetchNext(parfeval_futures, 3000);
        
        % Add results to list of results
        transfer(job_idx+1) = transfer_job;

        % Update progress plot
        if i == 1
            ax = mod_prog(epsilon, transfer);
        else
            ax = mod_prog(epsilon, transfer, ax);
        end
        
        % If a canidate flip is found
        if any(transfer > 0) && any(transfer < 0)
            epsilon(transfer == 0) = [];
            transfer(transfer == 0) = [];

            % Sort epsilon into correct order
            [epsilon, idx] = sort(epsilon);
            transfer = transfer(idx);
            
            % Search for neighboring pair of flip
            for j = 1:length(epsilon)-1
                if transfer(j)*transfer(j+1) < 0 
                    
                    % Cancel any currently running or queued jobs
                    unfinished = arrayfun(@(x) strcmp(x.State,'queued') || ...
                        strcmp(x.State,'running'),  parfeval_futures);
                    parfeval_futures = parfeval_futures(unfinished);
                    cancel(parfeval_futures);
                    
                    % Retest values at canidate a much higher accuracy 
                    test_futures(1) = parfeval(@optimal_rotation, 1, epsilon(j), ...
                        C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, 64000);
                    
                    test_futures(2) = parfeval(@optimal_rotation, 1, epsilon(j+1), ...
                        C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, 64000);
                    
                    % Retreive result
                    [~, test_job(1)] = fetchNext(test_futures);
                    [~, test_job(2)] = fetchNext(test_futures);
                    
                    % If flip still confirmed exit with results
                    if test_job(1)*test_job(2) < 0 
                        found_flip = true;
                        epsilon_low = epsilon(j);
                        epsilon_high = epsilon(j+1);
                        return;
                    else
                        if isempty(parfeval_futures)
                            epsilon_low = epsilon(j);
                            epsilon_high = epsilon(j+1);
                            return;
                        end
                        % If not restart jobs
                        for k = 1:length(parfeval_futures)
                            parfeval_futures(k) = parfeval(parfeval_futures(k).Function, ...
                                parfeval_futures(k).NumOutputArguments, ...
                                parfeval_futures(k).InputArguments{:});
                        end
                    end
                end
            end
        end
    end
else
    for i = 1:line_range      
        % Test points extending out from initial guess
        if mod(i,2) 
            flip = 1;
        else
            flip = -1;
        end
        amp = floor(i+1/2)/128;
        epsilon(i+1) = epsilon_0*(1-amp*flip);
        transfer(i+1) = optimal_rotation(epsilon(i+1), C, L, Q, num_modes, ...
            basis_modes, lambda, modal_amp, t_scale, tspan, init, 6000);
        
        % Update progress plot
        if i == 1
            ax = mod_prog(epsilon, transfer);
        else
            ax = mod_prog(epsilon, transfer, ax);
        end
        
        % If a canidate flip is found
        if any(transfer > 0) && any(transfer < 0)
            epsilon(transfer == 0) = [];
            transfer(transfer == 0) = [];

            % Sort into correct order
            [epsilon, idx] = sort(epsilon);
            transfer = transfer(idx);
            for j = 1:length(epsilon)-1
                if transfer(j)*transfer(j+1) < 0 
                    
                    % Retest values at canidate a much higher accuracy
                    test(1) = optimal_rotation(epsilon(j+1), C, L, Q, num_modes, ...
                         basis_modes, lambda, modal_amp, t_scale, tspan, init, 6000);
                    test(2) = optimal_rotation(epsilon(j+1), C, L, Q, num_modes, ...
                         basis_modes, lambda, modal_amp, t_scale, tspan, init, 6000);
                     
                     % If flip still confirmed exit with results
                    if test(1)*test(2) < 0 
                        found_flip = true;
                        epsilon_low = epsilon(j);
                        epsilon_high = epsilon(j+1);
                        return;
                    end
                end
            end
        end
    end
end

% Deal with edge cases
[epsilon, idx] = sort(epsilon);
transfer = transfer(idx);
[~, idx] = min(abs(transfer));
if idx > 2 && idx < length(epsilon)-1
    epsilon_low = epsilon(idx-1);
    epsilon_high = epsilon(idx+1);
elseif idx == 1
    epsilon_low = epsilon(1);
    epsilon_high = epsilon(2);
else
    epsilon_low = epsilon(end-1);
    epsilon_high = epsilon(end);
end
