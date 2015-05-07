function handle = plot_energy(modal_amp, t, direct, type, custom, id)
% Plot system energy from the time response of the Galerkin sytem
% Create figure handle and axis handle

h = figure;
ax = newplot;

num_modes = size(modal_amp,2);

TKE = sum(1/2*modal_amp'.^2);
plot(ax, t, TKE);
ax.XLabel.String = 'time (s)';
ax.YLabel.String = 'System Energy';
ax.Title.String = ['Predicted System Energy ' id];

legend(ax, ['System energy ' num2str(num_modes) ' modes']);

if numel(t) < 2 
   if nargout == 1
       handle = h;
   end
   return
end
Hz = 1/(t(2) - t(1));

if custom
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes) '_custom'];
else
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes)];
end

if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end
file_name = [direct_ext filesep 'Energy_' id '_'];

file_name = [file_name '_t' num2str(ceil(t(1))) '_' num2str(ceil(t(end))) ...
    's_' num2str(ceil(Hz)) 'Hz'];
drawnow;

% Save figure in Figure\Galerkin folder
saveas(h, file_name, 'fig');

if nargout == 1
    handle = h;
end
end