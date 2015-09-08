function futures = save_galerkin(direct, custom, time_int, calc_coef, classify_sim, ...
    i, futures, num_modes, results_coef, results_int, results_scores)
% SAVE_GALKERIN save all coefficients from the requested computation, saves all
% values on a separate thread so system isn't bottlenecked by harddrive
%
% futures = save_coef(direct, custom, time_int, calc_coef, i, futures,
% results_coef, results_int);

% If the last save hasn't finished wait
if i ~= 1
    wait(futures)
end

folder = 'Galerkin Coeff';

% Save all variables that were processed
pool = gcp('nocreate');

if ~isempty(pool)
    if time_int && calc_coef && classify_sim
        futures = parfeval(pool, @save_results, 0, num_modes, direct, folder, custom, ...
            results_coef, results_int, results_scores);
    elseif calc_coef && time_int
        futures = parfeval(pool, @save_results, 0, num_modes, direct, folder, custom, ...
            results_coef, results_int);
    elseif time_int && classify_sim
        futures = parfeval(pool, @save_results, 0, num_modes, direct, folder, custom, ...
            results_int, results_scores);
    elseif time_int 
        futures = parfeval(pool, @save_results, 0, num_modes, direct, folder, custom, results_int);
    else
        futures = parfeval(pool, @save_results, 0, num_modes, direct, folder, custom, results_coef);
    end
else
    if time_int && calc_coef && classify_sim
        save_results(num_modes, direct, folder, custom, ...
            results_coef, results_int, results_scores);
    elseif calc_coef && time_int
        save_results(num_modes, direct, folder, custom, ...
            results_coef, results_int);
    elseif time_int && classify_sim
        save_results(num_modes, direct, folder, custom, ...
            results_int, results_scores);
    elseif time_int 
        save_results(num_modes, direct, folder, custom, results_int);
    else
        save_results(num_modes, direct, folder, custom, results_coef);
    end
end
end