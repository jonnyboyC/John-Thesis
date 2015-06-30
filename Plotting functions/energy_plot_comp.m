function h = energy_plot_comp(results_pod, results_int, modes, models, legend_entries)

modal_amp = results_pod.modal_amp;
modal_amp_sim = results_int.modal_amp_sim;
t_scale = results_pod.u_scale/results_pod.l_scale;
t = results_int.t;

emp_TKE = sum(modal_amp(:,modes)'.^2)/2;
mean_emp_TKE = mean(emp_TKE);
std_emp_TKE = std(emp_TKE);
one_sec = [0 1];

h = figure;
h.Position = [400, 400, 650, 500];
ax = newplot;
hold(ax, 'on');

for i = 1:length(models)
    modal_amp_temp = modal_amp_sim{models(i), 1};
    t_temp = t{models(i), 1};%/t_scale;
    sim_TKE = sum(modal_amp_temp(:,2:end)'.^2)/2;
    plot(ax, t_temp, sim_TKE);
end

leg_names = cell(length(models)+2,1);

for i = 1:length(models)
    leg_names{i} = [num2str(length(modes)) ' modes ' legend_entries{i}];
end

leg_names{length(models)+1} = 'empirical mean';
leg_names{length(models)+2} = 'empirical 3 std';


plot(one_sec, [mean_emp_TKE, mean_emp_TKE], 'k');
plot(one_sec, 1*[std_emp_TKE, std_emp_TKE] + mean_emp_TKE, 'k-.');
plot(one_sec, -1*[std_emp_TKE, std_emp_TKE] + mean_emp_TKE, 'k-.');

legend(ax, leg_names, 'Location', 'BestOutside', 'Orientation', 'vertical');
ax.YScale = 'log';
ax.XLim = [0 1];
ax.YLim = [1/10000 100];
ax.YLabel = ylabel('System Energy', 'FontName', 'Time New Roman', 'FontSize', 14);
ax.XLabel = xlabel('Time(s)', 'FontName', 'Time New Roman', 'FontSize', 14);
ax.Title = title('Predicted Total Kinetic Energy', 'FontName', 'Time New Roman', 'FontSize', 14);
saveas(h, [num2str(length(modes)) '_' legend_entries{1}], 'png');