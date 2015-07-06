function futures = save_galerkin(direct, custom, time_int, calc_coef, classify_sim, ...
    i, futures, num_modes, results_coef, results_int, results_scores)
% SAVE_COEF save all coefficients from the request computation, saves all
% values on a separate thread so system isn't bottlenecked by harddrive
%
% futures = save_coef(direct, custom, time_int, calc_coef, i, futures,
% results_coef, results_int);

% If the last save hasn't finished wait
if i ~= 1
    wait(futures)
end

% Save all variables that were processed
pool = gcp;
if time_int && calc_coef && classify_sim
    futures = parfeval(pool, @save_results, 0, num_modes, direct, custom, results_coef, ...
        results_int, results_scores);
elseif calc_coef && time_int
    futures = parfeval(pool, @save_results, 0, num_modes, direct, custom, results_coef, ...
        results_int);
elseif time_int && classify_sim
    futures = parfeval(pool, @save_results, 0, num_modes, direct, custom, results_int, ...
        results_scores);
elseif time_int 
    futures = parfeval(pool, @save_results, 0, num_modes, direct, custom, results_int);
else
    futures = parfeval(pool, @save_results, 0, num_modes, direct, custom, results_coef);
end
end

function save_results(num_modes, direct, custom, varargin)
% check if folder exist create if empty
if custom
    direct_ext = [direct filesep 'Galerkin Coeff' filesep 'modes_' num_modes '_custom'];
else
    direct_ext = [direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(num_modes)];
end

if ~exist(direct_ext, 'dir')
    mkdir(direct_ext);
end

% file name for mat file
filename = [direct_ext filesep 'Coefficients_run_' num2str(varargin{1}.run_num) '.mat'];

% allow direct access to file
file = matfile(filename, 'Writable', true);

% save all 
for i = 1:length(varargin)
    file.(varargin{i}.name) = varargin{i};
end
end