function handle = plot_amp(modal_amp, t, direct, init)
% Create figure handle and axis handle
h = figure;
ax = newplot;
plot(ax, t, modal_amp);
ax.XLabel.String = ('time (s)');
ax.YLabel.String = ('Modal Amplitude');
ax.Title.String = ('Predicted Modal Amplitudes');

% Add amplitude legend
leg_names = cell(size(modal_amp, 2), 1);
for i = 1:size(modal_amp,2);
    leg_names{i} = ['Modal Amplitude ' num2str(i)]; 
end
legend(leg_names);
Hz = 1/(t(2) - t(1));

% Save figure in Figure\Galerkin folder
saveas(h, [direct '\Figures\Galerkin\Galerkin_' num2str(size(modal_amp,2)) '_modes_init' ...
    num2str(init) '_t' num2str(ceil(t(1))) '_' num2str(ceil(t(end))) 's_' ...
    num2str(ceil(Hz)) 'Hz'], 'fig');
drawnow;

if nargout == 1
    handle = h;
end

