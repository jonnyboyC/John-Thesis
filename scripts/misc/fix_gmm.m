function fix_gmm(varargin)

if nargin == 0
    direct = '';
else
    direct = varargin{1};
end

% Prompt the user for location of folder if not provided
fprintf(1, 'Please choose test data directory\n');
if strcmp(direct, '')
    direct = uigetdir('', 'Choose Source Image Directory');
end

try
    [direct_POD, ~] = prompt_folder('POD', 'first', direct);
catch
    fprintf('POD_Gen has not been run for %s\n', direct);
    return;
end

vars = matfile(direct_POD, 'Writable', true);

if isfield(vars, 'results_clust')
    results_clust = vars.results_clust;
    if isfield(results_clust, 'gm')
        results_clust.gmm = results_clust.gm;
        results_clust = rmfield(results_clust, 'gm');
        fprintf('Data directory %s fixed gm to gmm\n', direct);
    end
    vars.results_clust = results_clust;
else
    fprintf('Clustering has not been performed on %s%', direct);
end
