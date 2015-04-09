function plot_prediction(pod_u, pod_v, x, y, bnd_idx, modal_amp, t, dimensions, direct, id)
% Create a movie of the time response of the predicted Galerkin sytem
plot_points = 1:size(modal_amp,1);
modes = size(modal_amp,2);

% Fill all plots with blank images to set renderer to opengl
h_parent = figure('Name', ['Time Prediction for ' id], ...
           'Position', [500, 500, 500, 400]);

set(h_parent, 'Renderer', 'opengl');

if length(t) < 2 
    return;
end
Hz = 1/(t(2) - t(1));

% Intialize Video creator
if ~exist([direct filesep 'Figures' filesep 'Movies' filesep 'modes_' num2str(modes)], 'dir') 
    mkdir([direct filesep 'Figures' filesep 'Movies' filesep 'modes_' num2str(modes)]);
end

file_name = [direct filesep 'Figures' filesep 'Movies' filesep 'modes_' num2str(modes) filesep 'Flow_prediction_'];
if nargin == 12
    file_name = [file_name id];
end
writer = VideoWriter([file_name '_' num2str(t(1)) '_' num2str(t(end)) 's_' num2str(ceil(Hz)) 'Hz.avi']);
writer.Quality = 100;
writer.FrameRate = 60;
open(writer);

%% TODO May need to relook at this to make it more memory efficient
data_temp.xg = x;
data_temp.yg = y;
data_temp.pod = zeros(dimensions(1), dimensions(2));

data_u = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_v = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_m = zeros(dimensions(1), dimensions(2), size(plot_points,2));

% calculate predicted images
idx = 1;
sum_j = 1:modes;
for i = plot_points;
    data_u(:,:,idx) =  data_u(:,:,idx)...
        + reshape(pod_u(:,sum_j)*modal_amp(i,sum_j)',dimensions(1), dimensions(2));
    data_v(:,:,idx) =  data_v(:,:,idx)...
        + reshape(pod_v(:,sum_j)*modal_amp(i,sum_j)',dimensions(1), dimensions(2));
    data_m(:,:,idx) = sqrt(data_u(:,:,idx).^2+data_v(:,:,idx).^2);
    idx = idx + 1;
end

type = {'Flow Visualization'};

% Determine max values for u v and magnitude
% Leave a bit left off so useful plots can be produced even with blow up
[~, idx] = sort(abs(data_m(:)));
cmax = abs(data_m(idx(0.95*floor(length(idx)))));
cmin = -cmax;

data_temp.x = x;
data_temp.y = y;
data_temp.bnd_idx = bnd_idx;

% Plot results, print current image number, and save images to .avi video
for i = 1:size(plot_points,2)  
    fprintf('image %d of %d\n', i, size(plot_points,2));
    for j = 1:length(type);
        data_temp.pod = squeeze(data_m(:,:,i));
        data_temp.u = squeeze(data_u(:,:,i));
        data_temp.v = squeeze(data_v(:,:,i));
        if i == 1
            figure(h_parent);
            [h_child, ax] = Plottec2(data_temp);
            ax.Title = title(type{j}, 'fontname','times new roman','fontsize', 14);
            ax.XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
            ax.YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
            ax.ZLim = [cmin, cmax];
            ax.CLim = [cmin, cmax];
            colorbar;
        else
            h_child = Plottec2(data_temp, h_child);
        end
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end