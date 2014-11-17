function h = amp_comp(modal_amp, modal_amp_til, t0, eig_func, tspan, mode)
close all;

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
t = t0(t0 < tspan(2) & t0 > tspan(1));
modal_amp = modal_amp((t0 < tspan(2) & t0 > tspan(1)), mode);
plot(ax, t, modal_amp, 'r');

modal_amp_til = modal_amp_til((t0 < tspan(2) & t0 > tspan(1)), mode);
plot(ax, t, modal_amp_til, 'g');

ax.Title.String = ['Predicted Modal Amplitude a_' num2str(mode)];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Amplitude';
hold off

end