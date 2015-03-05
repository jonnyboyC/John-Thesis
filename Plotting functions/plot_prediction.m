function plot_prediction(pod_u, pod_v, pod_vor, x, y, bnd_idx, modal_amp, t, num_pods, dimensions, direct, id)
% Create a movie of the time response of the predicted Galerkin sytem

% Check to match sure requested image isn't too large
% if size(modal_amp,1) > 5000
%     step_size = ceil(size(modal_amp,1)/5000);
%     plot_points = 1:step_size:size(modal_amp,1);
% else
    plot_points = 1:size(modal_amp,1);
% end

% Fill all plots with blank images to set renderer to opengl
dummie = zeros(2,2);
h = figure('Name', ['Time Prediction for ' id], ...
           'Position', [500, 500, 800, 400]);
subplot(2,2,1)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

subplot(2,2,2)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

subplot(2,2,3)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

subplot(2,2,4)
pcolor(dummie);
axis tight
set(gca, 'nextplot', 'replacechildren');
set(h, 'Renderer', 'opengl');

if length(t) < 2 
    return;
else
Hz = 1/(t(2) - t(1));

% Intialize Video creator
ext = '\Figures\Movies\POD_';
if nargin == 10
    ext = [ext id];
end
writer = VideoWriter([direct ext num2str(num_pods) '_' ...
   num2str(t(1)) '_' num2str(t(end)) 's_' num2str(ceil(Hz)) 'Hz_Galerkin.avi']);
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
data_vor = zeros(dimensions(1), dimensions(2), size(plot_points,2));

% calculate predicted images
idx = 1;
for i = plot_points;
    for j = 1:num_pods
        data_u(:,:,idx) =  data_u(:,:,idx)...
            + reshape(pod_u(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
        data_v(:,:,idx) =  data_v(:,:,idx)...
            + reshape(pod_v(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
        data_vor(:,:,idx) =  data_vor(:,:,idx)...
            + reshape(pod_vor(:,j)*modal_amp(i,j),dimensions(1), dimensions(2));
    end
    data_m(:,:,idx) = sqrt(data_u(:,:,idx).^2+data_v(:,:,idx).^2);
    idx = idx + 1;
end

data = {data_u, data_v, data_m, data_vor};
type = {'Streamwise Velocity', 'Spanwise Velocity', 'Velocity Magnitude', 'Vorticity'};

% Determine max values for u v and magnitude
cmax = zeros(size(data,2),1);
cmin = zeros(size(data,2),1);
for i = 1:length(data)
    cmax(i) = max(abs(data{i}(:)));
    cmin(i) = -cmax(i);
end

% Preallocated figure handles and axes handles
h_sub = gobjects(4,1);
ax_sub = gobjects(4,1);

data_temp.x = x;
data_temp.y = y;

% Plot results, print current image number, and save images to .avi video
for i = 1:size(plot_points,2)  
    fprintf('image %d of %d\n', i, size(plot_points,2));
    for j = 1:length(data);
        data_temp.pod = squeeze(data{j}(:,:,i));
        data_temp.cmax = cmax(j);
        subplot(2,2,j);
        if i == 1
            [h_sub(j), ax_sub(j)] = Plottec2(data_temp, 0, bnd_idx);
            ax_sub(j).Title = title(type{j}, 'fontname','times new roman','fontsize', 14);
            ax_sub(j).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
            ax_sub(j).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
            colorbar;
        else
            h_sub(j) = Plottec2(data_temp, h_sub(j));
        end
        ax_sub(j).ZLim = [cmin(j), cmax(j)];
        ax_sub(j).CLim = [cmin(j), cmax(j)];
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end