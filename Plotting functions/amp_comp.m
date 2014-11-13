function h = amp_comp(modal_amp, t, modal_amp_til, t_til, eig_func, tspan, mode)

mode_dev1 = std(eig_func(:,mode));
mode_dev2 = mode_dev1*2;
mode_mean = mean(eig_func(:,mode));

h = figure;
ax = newplot;

plot(ax, tspan, [mode_dev1 + mode_mean, mode_dev1 + mode_mean], 'b-.', ...
         tspan, [-mode_dev1 + mode_mean, -mode_dev1 + mode_mean] ,'b-.');
hold on
plot(ax, tspan, [mode_dev2 + mode_mean, mode_dev2 + mode_mean], 'k-.', ...
         tspan, [-mode_dev2 + mode_mean, -mode_dev2 + mode_mean] ,'k-.');
t = t(t < tspan(2));
modal_amp = modal_amp((t < tspan(2)), mode);
plot(ax, t, modal_amp, 'r');

t_til = t(t_til < tspan(2));
modal_amp_til = modal_amp_til((t_til < tspan(2)), mode);
plot(ax, t_til, modal_amp_til, 'g');

ax.Title.String = ['Predicted Modal Amplitude a_' num2str(mode)];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude';
hold off

end