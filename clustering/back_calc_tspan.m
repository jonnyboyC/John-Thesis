function [tspan, multiplier] = back_calc_tspan(exp_sampling_rate, integration, modal_amp)
% BACK_CALC_TSPAN because I am dumb and hadn't been storing this before
% have to not back calculate the necessary values, should be able to
% completely eliminate this in the future

m = flow_comps(integration.t);
models = flow_ncomps(integration.t);
max_steps = size(modal_amp,1);
multiplier = 0;

for i = 1:models;
    s = flow_comps(integration.t.(m{i}));
    sub_models = flow_ncomps(integration.t.(m{i}));
    
    for j = 1:sub_models
        time_steps = length(integration.t.(m{i}).(s{j}));
        if  multiplier == 0 && time_steps > 1
            t = integration.t.(m{i}).(s{j});
            multiplier = 1/((t(2)-t(1))*exp_sampling_rate);
        end
        if time_steps >= max_steps
            if time_steps == max_steps
                tspan = integration.t.(m{i}).(s{j});
                return;
            end
            if time_steps > max_steps
                max_steps = time_steps;
            end
        end
    end
end

% if tspan not assigned it is only used for a message
tspan = 0;

