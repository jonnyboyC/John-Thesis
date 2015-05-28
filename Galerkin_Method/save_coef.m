function futures = save_coef(direct, custom, time_int, calc_coef, i, futures, results_coef, results_int)
% SAVE_COEF save all coefficients from the request computation, saves all
% values on a separate thread so system isn't bottlenecked by harddrive
%
% futures = save_coef(direct, custom, time_int, calc_coef, i, futures,
% results_coef, results_int);

if i ~= 1
    wait(futures)
end
pool = gcp;
if time_int && calc_coef
    futures = parfeval(pool, @save_results, 0, direct, true, custom, results_coef, ...
        results_int);
elseif time_int
    futures = parfeval(pool, @save_results, 0, direct, false, custom, results_int);
else
    futures = parfeval(pool, @save_results, 0, direct, true, custom, results_coef);
end
end

function save_results(direct, coef, custom, varargin)
% check if folder exist create if empty
if custom
    direct_ext = [direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(varargin{1}.num_modesG) '_custom'];
else
    direct_ext = [direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(varargin{1}.num_modesG)];
end

if ~exist(direct_ext, 'dir')
    mkdir(direct_ext);
end

% file name for mat file
filename = [direct_ext filesep 'Coefficients_run_' num2str(varargin{1}.run_num) '.mat'];

% allow direct access to file
file = matfile(filename, 'Writable', true);

% if only results_coef was passed only save those
if nargin == 4
    if coef
        file.results_coef = varargin{1};
    else
        file.results_int = varargin{1};
    end
else
    file.results_coef = varargin{1};
    file.results_int = varargin{2};
end
end