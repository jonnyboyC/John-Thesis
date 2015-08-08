modes = [2:3, 5:17];
models = {[1,2], [6, 8, 12, 14], [5, 7, 11, 13], [2, 3, 9 10]};
legend = {{'GM', 'GM_W'}, ...
          {'GM1', 'GM1_W', 'GM1_N', 'GM1_{NW}'}, ...
          {'GM2', 'GM2_W', 'GM2_N', 'GM2_{NW}'}, ...
          {'GM3', 'GM3_W', 'GM3_N', 'GM3_{NW}'}};
for i = 1:4
    energy_plot_comp(results_pod, results_int, modes, models{i}, legend{i});
end