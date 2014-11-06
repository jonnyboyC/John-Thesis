function handle = plot_amp(modal_amp, t, h)
ax = newplot;
plot(ax, t, modal_amp);
ax.XLabel.String = ('time (s)');
ax.YLabel.String = ('Modal Amplitude');
ax.Title.String = ('Predicted Modal Amplitudes');

leg_names = cell(size(modal_amp, 2), 1);
for i = 1:size(modal_amp,2);
   leg_names{i} = ['Modal Amplitude ' num2str(i)']; 
end
legend(leg_names);
drawnow;

if nargout == 1
    handle = h;
end

