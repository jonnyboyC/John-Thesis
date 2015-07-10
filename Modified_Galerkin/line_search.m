function [epsilon_low, epsilon_high, transfer, flip] = line_search(problem)
% LINE_SEARCH brute force sweep of an area around the initial condition
% suggested by balajawecz.
%
% [epsilon_low, epsilon_high, transfer, flip] = LINE_SEARCH(problem) edit
% line_search for input details

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

% Find initial value for line search and plot
transfer(1) = optimal_rotation(epsilon_0, C, L, Q, num_modes, basis_modes, lambda, modal_amp, t_scale, tspan, init, 6000);

% Brute force Line Search, rough pass parrallel loop. Idea is to quickly
% sweep a large area to look for sign changes to use for a finer pass in
% fzero


% prepare asynounous jobs to run since the number of iterations may vary
% notably
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
    [job_idx, transfer_job] = fetchNext(parfeval_futures, 3000);
    transfer(job_idx+1) = transfer_job;
    
    if i == 1
        ax = mod_prog(epsilon, transfer);
    else
        ax = mod_prog(epsilon, transfer, ax);
    end
    if any(transfer > 0) && any(transfer < 0)
        epsilon(transfer == 0) = [];
        transfer(transfer == 0) = [];
        
        [epsilon, idx] = sort(epsilon);
        transfer = transfer(idx);
        for j = 1:length(epsilon)-1
            if transfer(j)*transfer(j+1) < 0 
                flip = true;
                epsilon_low = epsilon(j);
                epsilon_high = epsilon(j+1);
                cancel(parfeval_futures)
                return;
            end
        end
    end
end

flip = false;

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
