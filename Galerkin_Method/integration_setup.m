function system = integration_setup(system, i, num_modes)
% INTEGRATION_SETUP convert the system coefficients from thier logical form
% to that used by the system solve TIME_INTEGRATION
%
%   [reduced_model_coeff] = integration_setup(eddy, vis, l, q, i,
%       total_models, linear_models, num_modes);

Mod = false;

s = flow_comps(system{i}.l);
s_comps = flow_ncomps(system{i}.l);

e = flow_comps(system{i}.eddy);
e_comps = flow_ncomps(system{i}.eddy);

vis = system{i}.vis;

linear_nonlinear = 1:2;

% Generate the coef structure, add in spots for nonlinear 
for j = linear_nonlinear
    for k = 1:e_comps
        for kk = 1:s_comps
            
            % Pre-combine eddy visocity if we are using a linear system
            if j == 1
                total_vis = repmat(system{i}.eddy.(e{k}).(s{kk}) + vis, 1, num_modes);
                coeff = [system{i}.l.(s{kk}).*total_vis, system{i}.q.(s{kk})];
            else
                if k == 1
                    continue;
                end
                coeff = [system{i}.l.(s{kk}), system{i}.q.(s{kk})];
            end
            
            % Get a name for each model
            if j == 1
                type = s{kk};
            else
                type = ['NL_' s{kk}];
                system{i}.eddy.(e{k}).(type) = system{i}.eddy.(e{k}).(s{kk});
            end

            % Move all coefficients into one row for each mode
            system{i}.coef.(e{k}).(type) = ode_coefficients(num_modes, coeff, Mod);
        end
    end
end
    
end