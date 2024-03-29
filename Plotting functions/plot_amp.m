function handle = plot_amp(modal_amp, t, direct, type, custom, id)
% Plot modal amplitudes from the time response of the Galerkin sytem
% Create figure handle and axis handle
h = figure;
ax = newplot;

num_modes = size(modal_amp,2);

plot(ax, t, modal_amp);
ax.XLabel.String = 'time (s)';
ax.YLabel.String = 'Modal Amplitude';
ax.Title.String = ['Predicted Modal Amplitudes ' id];
axis(ax, 'tight');


% Add amplitude legend
leg_names = cell(size(modal_amp, 2), 1);
for i = 1:size(modal_amp,2);
    leg_names{i} = ['Modal Amplitude ' num2str(i)]; 
end
legend(ax, leg_names);
legend(ax, 'hide')
if numel(t) < 2 
   if nargout == 1
       handle = h;
   end
   return
end
Hz = 1/(t(2) - t(1));

if custom
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes) '_custom' filesep 'amplitude'];
else
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes) filesep 'amplitude'];
end

if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end
file_name = [direct_ext filesep 'Amplitude_' strrep(id, ' ', '_')];

file_name = [file_name '_t' num2str(ceil(t(1))) '_' num2str(ceil(t(end))) ...
    's_' num2str(ceil(Hz)) 'Hz'];
drawnow;

% Save figure in Figure\Galerkin folder
saveas(h, file_name, 'fig');
saveas(h, file_name, 'png');

if nargout == 1
    handle = h;
end
end

