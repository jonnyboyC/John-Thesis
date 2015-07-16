function system_coefs = load_coef(direct, run_num, num_modes, custom)
% LOAD_COEF load the coefficients of a previous run for time integration
%
% system_coefs = LOAD_COEF(direct, run_num, custom)

% load from custom folder if requested
if custom
    vars = load([direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(num_modes-1) '_custom' ...
        filesep 'Coefficients_run_' num2str(run_num) '.mat'], 'results_coef');
else
    vars = load([direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(num_modes-1) ...
        filesep 'Coefficients_run_' num2str(run_num) '.mat'], 'results_coef');
end

system_coefs = vars.results_coef.system;
end