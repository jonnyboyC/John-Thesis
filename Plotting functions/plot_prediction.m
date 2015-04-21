function plot_prediction(pod_u, pod_v, x, y, bnd_idx, modal_amp, t, dimensions, direct, id)
% Create a movie of the time response of the predicted Galerkin sytem
max_plot = 500;
if max_plot < size(modal_amp,1);  
    plot_points = 1:max_plot;
else
    plot_points = 1:size(modal_amp,1);
end
modes = size(modal_amp,2);

% Fill all plots with blank images to set renderer to opengl
h_parent = figure('Name', ['Time Prediction for ' id], ...
           'Position', [500, 400, 700, 500]);

set(h_parent, 'Renderer', 'opengl');

if length(t) < 2 
    return;
end
Hz = 1/(t(2) - t(1));

% Intialize Video creator
if ~exist([direct filesep 'Figures' filesep 'Movies' filesep 'modes_' num2str(modes-1)], 'dir') 
    mkdir([direct filesep 'Figures' filesep 'Movies' filesep 'modes_' num2str(modes-1)]);
end

file_name = [direct filesep 'Figures' filesep 'Movies' filesep 'modes_' ...
             num2str(modes-1) filesep 'Flow_prediction_' id];

% TODO May need to relook at this to make it more memory efficient
data_u = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_v = zeros(dimensions(1), dimensions(2), size(plot_points,2));

% calculate predicted images
sum_j = 1:modes;
for i = plot_points;
    data_u(:,:,i) = reshape(pod_u(:,sum_j)*modal_amp(i,sum_j)',dimensions(1), dimensions(2));
    data_v(:,:,i) = reshape(pod_v(:,sum_j)*modal_amp(i,sum_j)',dimensions(1), dimensions(2));
end

data_m = sqrt(data_u.^2 + data_v.^2);


type = {'Flow Visualization'};

% Determine max values for u v and magnitude
% Leave a bit left off so useful plots can be produced even with blow up
[~, idx] = sort(abs(data_m(:)));
cmax = abs(data_m(idx(floor(0.95*length(idx)))));
cmin = 0;

data_temp.x = x;
data_temp.y = y;
data_temp.bnd_idx = bnd_idx;

writer = VideoWriter([file_name '_' num2str(t(1)) '_' num2str(t(end)) 's_' num2str(ceil(Hz)) 'Hz.avi']);
writer.Quality = 100;
writer.FrameRate = 60;
open(writer);

% Plot results, print current image number, and save images to .avi video
for i = 1:size(plot_points,2)  
    fprintf('image %d of %d\n', i, size(plot_points,2));
    for j = 1:length(type);
        data_temp.pod = squeeze(data_m(:,:,i));
        data_temp.u = squeeze(data_u(:,:,i));
        data_temp.v = squeeze(data_v(:,:,i));
        if i == 1
            figure(h_parent);
            [h_surf, h_quiver, ax] = plot_flow(data_temp);
            ax.Title = title(type{j}, 'fontname','times new roman','fontsize', 14);
            ax.XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
            ax.YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
            ax.ZLim = [cmin, cmax];
            ax.CLim = [cmin, cmax];
            colorbar;
        else
            [h_surf, h_quiver] = plot_flow(data_temp, h_surf, h_quiver);
        end
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end