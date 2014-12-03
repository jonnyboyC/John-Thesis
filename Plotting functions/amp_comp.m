function h = amp_comp(coeff, tspan, eig_func, num_pods, mode2plot, init)
% Plot the time response of Galerkin System verse the std dev of the direct
% PIV projection onto the POD basis

% Calculate standard deviation and mean
mode_dev1 = std(eig_func(:,mode2plot))*3;
mode_dev2 = mode_dev1*2;
mode_mean = mean(eig_func(:,mode2plot));

h = figure;
ax = newplot;

% Plot first and second standard deviation line
plot(ax, [tspan(1) tspan(end)], [mode_dev1 + mode_mean, mode_dev1 + mode_mean], 'b-.', ...
         [tspan(1) tspan(end)], [-mode_dev1 + mode_mean, -mode_dev1 + mode_mean] ,'b-.');
hold on
plot(ax, [tspan(1) tspan(end)], [mode_dev2 + mode_mean, mode_dev2 + mode_mean], 'k-.', ...
         [tspan(1) tspan(end)], [-mode_dev2 + mode_mean, -mode_dev2 + mode_mean] ,'k-.');
     
% Solve for common time scale
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
colors = {'b-', 'c-', 'g-'};
for i = 1:size(coeff,2)
	reduced_model_coeff = -ode_coefficients(num_pods, num_pods, coeff{i});
    tic1 = tic;
    [t, modal_amp] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff), tspan, ...
        eig_func(init,1:num_pods), options);
    toc(tic1);
    plot(ax, t, modal_amp(:,mode2plot), colors{i});
end
% legend(ax, 'ROM1', 'ROM2', 'ROM3');
ax.Title.String = ['Predicted Modal Amplitude a_' num2str(mode2plot)];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude';
hold off

end