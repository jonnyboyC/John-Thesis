function [epsilon, transfer, flip_idx] = line_search(problem)

C           = problem.C;
L           = problem.L;
Q           = problem.Q;
OG_nm       = problem.OG_nm;
RD_nm       = problem.RD_nm;
lambda      = problem.lambda;
modal_amp   = problem.modal_amp;
t           = problem.t;
init        = problem.init;
line_range  = problem.line_range;
flip_idx    = 0;

% intial guess for epsilon
epsilon_0 = sum(L(1:RD_nm, 1:RD_nm)*lambda(1:RD_nm)); 
epsilon = epsilon_0;

% Find initial value for line search and plot
transfer = optimal_rotation(epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_amp, t, init, 18000);

% Brute force Line Search, rough pass parrallel loop. Idea is to quickly
% sweep a large area to look for sign changes to use for a finer pass in
% fzero

search_space = 1:20:line_range;
exclusion_list = [];
for i = 1:length(search_space)-1
    range        = search_space(i+1) - search_space(i);
    search_range = search_space(i):search_space(i+1)-1;
    
    epsilon_temp  = zeros(range, 1);
    transfer_temp = zeros(range, 1);
    parfor j = 1:20
        search = search_range(j)
        if mod(j, 2)
            flip = 1;
        else 
            flip = -1;
        end
        amp = floor((search)/2)/8;
        epsilon_temp(j) = epsilon_0*(1-amp*flip);
        transfer_temp(j) = optimal_rotation(epsilon_temp(j), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, t, init, 18000);
    end

    epsilon       = [epsilon; epsilon_temp];
    transfer      = [transfer; transfer_temp];

    % Get sorted list of epsilon and unresolved transfer terms to try to find
    % and interval of sign change
    [epsilon, idx] = sort(epsilon);
    transfer = transfer(idx);
    sign_t = sign(transfer);
    
    % Update search progress
    if i == 1
        ax = mod_prog(epsilon, transfer);
    else
        ax = mod_prog(epsilon, transfer, ax);
    end


    % Find suspect zeros perform higher interation rotation if potential
    % bound found
    for k = 1:length(sign_t)-1
        if sign_t(k)*sign_t(k+1) < 0 && ~any(exclusion_list == k)
            transfer1 = optimal_rotation(epsilon(k), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, t, init, 18000);
            transfer2 = optimal_rotation(epsilon(k+1), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, t, init, 18000);
            if transfer1*transfer2 < 0 
                flip_idx = k;
                return;
            else
                exclusion_list = [exclusion_list, k];
            end
        end
    end
end
end