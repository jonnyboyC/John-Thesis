function plot_prediction(pod_U, X, bnd_idx, modal_amp, t, dimensions, streamlines, direct, custom, id)
% Create a movie of the time response of the predicted Galerkin sytem
max_plot = 2000;
if max_plot < size(modal_amp,1);  
    plot_points = 1:max_plot;
else
    plot_points = 1:size(modal_amp,1);
end
num_modes = size(modal_amp,2);

% Fill all plots with blank images to set renderer to opengl
h_parent = figure('Name', ['Time Prediction for ' id], ...
           'Position', [0, 0, 700, 800]);
movegui(h_parent,'center')
set(h_parent, 'Renderer', 'opengl');

if length(t) < 2 
    return;
end
Hz = 1/(t(2) - t(1));

% Determine the file the location that the movie should be placed
if custom
    direct_ext = [direct filesep 'Figures' filesep 'Movies' filesep 'modes_' ...
        num2str(num_modes-1) '_custom'];
else
    direct_ext = [direct filesep 'Figures' filesep 'Movies' filesep 'modes_' ...
        num2str(num_modes-1)];
end

% Intialize Video creator
if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end

% 
file_name = [direct_ext filesep 'Flow_prediction_' strrep(id, ' ', '_')];

dims = flow_dims(X);
[~, u] = flow_comps_ip(X, pod_U);

% Shorts handles for modes to sum
sum_full = 1:num_modes;
sum_flux = 2:num_modes;

type = {'Full Flow Visualization', 'Turbulent Flow Visualization'};

data_temp.X = X;
data_temp.bnd_idx = bnd_idx;

writer = VideoWriter([file_name '_' num2str(t(1)) '_' num2str(t(end)) 's_' ...
    num2str(ceil(Hz)) 'Hz'], 'MPEG-4');
writer.Quality = 100;
writer.FrameRate = 60;
open(writer);

% Preallocate figure subplot handles
h_surf = gobjects(2, 1);
ax = gobjects(2, 1);
if streamlines 
    h_dir = cell(2, 1);
else
    h_dir = gobjects(2, 1);
end


% Plot results, print current image number, and save images to .avi video
for i = 1:size(plot_points,2)  
    fprintf('image %d of %d\n', i, size(plot_points,2));
    for j = 1:2
        if j == 1
            for k = 1:dims
                data_temp.U.(u{k}) = reshape(pod_U.(u{j})(:,sum_full)*modal_amp(i,sum_full)',dimensions);
            end
        else
            for k = 1:dims
                data_temp.U.(u{k}) = reshape(pod_U.(u{j})(:,sum_flux)*modal_amp(i,sum_flux)',dimensions);
            end
        end
        figure(h_parent);
        subplot(2,1,j);
        if i == 1
            [h_surf(j), h_dir(j), ax(j)] = plot_vector_field(data_temp, streamlines);
            ax(j).Title = title(type{j}, 'fontname','times new roman','fontsize',14);
            ax(j).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
            ax(j).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
            colorbar;
        else
            [h_surf(j), h_dir(j)] = plot_vector_field(data_temp, streamlines, h_surf(j), h_dir(j));
        end
    end
    frame = getframe(gcf);
    writeVideo(writer, frame);
end
end