function plot_prediction(pod_u, pod_v, pod_vor, x, y, bnd_idx, modal_amp, t, dimensions, direct, custom, id)
% Create a movie of the time response of the predicted Galerkin sytem
max_plot = 500;
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

file_name = [direct_ext filesep 'Flow_prediction_' id];

% TODO May need to relook at this to make it more memory efficient
data_u_full = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_v_full = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_vor_full = zeros(dimensions(1), dimensions(2), size(plot_points,2));

data_u_flux = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_v_flux = zeros(dimensions(1), dimensions(2), size(plot_points,2));
data_vor_flux = zeros(dimensions(1), dimensions(2), size(plot_points,2));

% calculate predicted images
sum_full = 1:num_modes;
sum_flux = 2:num_modes;
for i = plot_points;
    data_u_full(:,:,i) = reshape(pod_u(:,sum_full)*modal_amp(i,sum_full)',dimensions);
    data_v_full(:,:,i) = reshape(pod_v(:,sum_full)*modal_amp(i,sum_full)',dimensions);
    data_vor_full(:,:,i) = reshape(pod_vor(:,sum_full)*modal_amp(i,sum_full)',dimensions);

    
    data_u_flux(:,:,i) = reshape(pod_u(:,sum_flux)*modal_amp(i,sum_flux)',dimensions);
    data_v_flux(:,:,i) = reshape(pod_v(:,sum_flux)*modal_amp(i,sum_flux)',dimensions);
    data_vor_flux(:,:,i) = reshape(pod_vor(:,sum_flux)*modal_amp(i,sum_flux)',dimensions);
end

% data_m_full = sqrt(data_u_full.^2 + data_v_full.^2);
% data_m_flux = sqrt(data_u_flux.^2 + data_v_flux.^2);

type = {'Full Flow Visualization', 'Turbulent Flow Visualization'};

% Determine max values for u v and magnitude
% Leave a bit left off so useful plots can be produced even with blow up
[~, idx] = sort(abs(data_vor_full(:)));
cmax(1) = abs(data_vor_full(idx(floor(0.95*length(idx)))));
cmin(1) = -cmax(1);

[~, idx] = sort(abs(data_vor_flux(:)));
cmax(2) = abs(data_vor_flux(idx(floor(0.95*length(idx)))));
cmin(2) = -cmax(2);

data_temp.x = x;
data_temp.y = y;
data_temp.bnd_idx = bnd_idx;

writer = VideoWriter([file_name '_' num2str(t(1)) '_' num2str(t(end)) 's_' ...
    num2str(ceil(Hz)) 'Hz.avi'], 'MPEG-4');
writer.Quality = 100;
writer.FrameRate = 60;
open(writer);

% Preallocate figure subplot handles
h_surf = gobjects(2, 1);
h_quiver = gobjects(2, 1);
if pretty 
    h_stream = gobjects(2, 1);
end
ax = gobjects(2, 1);

% Plot results, print current image number, and save images to .avi video
for i = 1:size(plot_points,2)  
    fprintf('image %d of %d\n', i, size(plot_points,2));
    for j = 1:2
        if j == 1
            data_temp.pod = squeeze(data_vor_full(:,:,i));
            data_temp.u = squeeze(data_u_full(:,:,i));
            data_temp.v = squeeze(data_v_full(:,:,i));
        else
            data_temp.pod = squeeze(data_vor_flux(:,:,i));
            data_temp.u = squeeze(data_u_flux(:,:,i));
            data_temp.v = squeeze(data_v_flux(:,:,i));
        end
        figure(h_parent);
        subplot(2,1,j);
        if i == 1
            [h_surf(j), h_quiver(j), ax(j)] = plot_vector_field(data_temp);
            if pretty
                h_stream(j) = streamslice(data_temp.x', data_temp.y', data_temp.u', data_temp.v');
            end
            ax(j).Title = title(type{j}, 'fontname','times new roman','fontsize',14);
            ax(j).XLabel = xlabel('x/D', 'fontname','times new roman','fontsize',12);
            ax(j).YLabel = ylabel('y/D', 'fontname','times new roman','fontsize',12);
            ax(j).ZLim = [cmin(j), cmax(j)];
            ax(j).CLim = [cmin(j), cmax(j)];
            colorbar;
        else
            [h_surf(j), h_quiver(j)] = plot_vector_field(data_temp, h_surf(j), h_quiver(j));
            if pretty
                delete(h_stream(j))
                h_stream(j) = streamslice(data_temp.x', data_temp.y', data_temp.u', data_temp.v');
            end
        end
    end
    frame = getframe(gcf);
    writeVideo(writer, frame)
end
end