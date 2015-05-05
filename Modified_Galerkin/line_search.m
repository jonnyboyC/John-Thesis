function [epsilon_low, epsilon_high, transfer, flip] = line_search(problem)

C           = problem.C;
L           = problem.L;
Q           = problem.Q;
OG_nm       = problem.OG_nm;
RD_nm       = problem.RD_nm;
lambda      = problem.lambda;
modal_amp   = problem.modal_amp;
tspan       = problem.tspan;
init        = problem.init;
line_range  = problem.line_range;


% intial guess for epsilon
epsilon_0 = sum(L(1:RD_nm, 1:RD_nm)*lambda(1:RD_nm)); 
epsilon = zeros(line_range+1,1);
transfer = zeros(line_range+1,1);
epsilon(1) = epsilon_0;

% Find initial value for line search and plot
transfer(1) = optimal_rotation(epsilon_0, C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, 64000);

% Brute force Line Search, rough pass parrallel loop. Idea is to quickly
% sweep a large area to look for sign changes to use for a finer pass in
% fzero

parfeval_futures = parallel.FevalOnAllFuture;
for i = 1:line_range
    if mod(i,2) 
        flip = 1;
    else
        flip = -1;
    end
    amp = floor(i+1/2)/16;
    epsilon(i+1) = epsilon_0*(1-amp*flip);
    parfeval_futures(i) = parfeval(@optimal_rotation, 1, epsilon(i+1), ...
        C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, 64000);
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

% 
% search_space = 1:20:line_range;
% exclusion_list = [];
% for i = 1:length(search_space)-1
%     range        = search_space(i+1) - search_space(i);
%     search_range = search_space(i):search_space(i+1)-1;
%     
%     epsilon_temp  = zeros(range, 1);
%     transfer_temp = zeros(range, 1);
%     parfor j = 1:20
%         search = search_range(j)
%         if mod(j, 2)
%             flip = 1;
%         else 
%             flip = -1;
%         end
%         amp = floor((search)/2)/16;
%         epsilon_temp(j) = epsilon_0*(1-amp*flip);
%         transfer_temp(j) = optimal_rotation(epsilon_temp(j), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, 64000);
%     end
% 
%     epsilon       = [epsilon; epsilon_temp];
%     transfer      = [transfer; transfer_temp];
% 
%     % Get sorted list of epsilon and unresolved transfer terms to try to find
%     % and interval of sign change
%     [epsilon, idx] = sort(epsilon);
%     transfer = transfer(idx);
%     sign_t = sign(transfer);
%     
%     % Update search progress
%     if i == 1
%         ax = mod_prog(epsilon, transfer);
%     else
%         ax = mod_prog(epsilon, transfer, ax);
%     end
% 
% 
%     % Find suspect zeros perform more interations of optimal_rotation to more assuredly
%     % find bounds
%     for k = 1:length(sign_t)-1
%         if sign_t(k)*sign_t(k+1) < 0 && ~any(exclusion_list == k)
%             transfer1 = optimal_rotation(epsilon(k), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, 64000);
%             transfer2 = optimal_rotation(epsilon(k+1), C, L, Q, OG_nm, RD_nm, lambda, modal_amp, tspan, init, 64000);
%             % if higher accuracy produces bounded solutions exit and
%             % indicate flip location
%             if transfer1*transfer2 < 0 
%                 flip_idx = k;
%                 return;
%             % Otherwise add point to exclusion list
%             else
%                 exclusion_list = [exclusion_list, k];
%             end
%         end
%     end
% end
end