function [tspan, sample_freq, multiplier] = calc_tspan(tspan, exp_sampling_rate, modal_amp)

if length(tspan) == 3
   extender = tspan{3}; 
else
   extender = 1;
end

if length(tspan) >= 2
    multiplier = tspan{2};
else
    multiplier = 1;
end

steps = length(modal_amp) - 1;
sample_freq = exp_sampling_rate*multiplier;

tspan = 0:1/sample_freq:(((steps + multiplier)*extender)/sample_freq);