function [h_magnitude, cax] = plot_scalar_field(data, h_mag)
% PLOT_SCALAR_FIELD plot an individual scalar field with magnitude
% indicated by color
%
% [h_magnitude, cax] = PLOT_SCALAR_FIELD(data, h_mag)

% Create pseudo color plot, for flow visualization with included boundaries
if nargin == 1 || all(isfield(data, {'bnd_X'}))
    ax = newplot;    
    
    X = data.X;
    x = flow_comps_ns(X);
    dimensions = size(X.(x{1}));
    
    h_mag = surface(X.(x{1}), X.(x{2}), zeros(dimensions), data.field);
    h_mag.FaceColor = 'interp';
    h_mag.EdgeColor = 'none';
    
    if isfield(data, 'bnd_idx')
        % Get boundary index
        bnd_idx = data.bnd_idx;
        black = zeros([dimensions, 3]);

        % Plot boundary index
        hold(ax, 'on');
        h2 = surf(ax, X.(x{1}), X.(x{2}), zeros(dimensions), black);
        hold(ax, 'off');
        
        % Make all non-boundary location transparent
        bnd_idx = double((bnd_idx == -1));
        h2.AlphaDataMapping = 'none';
        h2.FaceColor = 'flat';
        h2.EdgeColor = 'none';
        h2.AlphaData = bnd_idx;
        h2.FaceAlpha = 'flat';
        
        if all(isfield(data, {'bnd_X'}))
            % Get Open Flow boundaries
            bnd_X = data.bnd_X;
            red = zeros([dimensions 3]);
            red(:,:,1) = ones(dimensions);
            
            % Plot open flow boundaries
            hold(ax, 'on');
            h3 = surf(ax, X.(x{1}), X.(x{2}), zeros(dimensions), red);
            hold(ax, 'off');
            
            % Make all non-boundary location transparent
            flow_boundary = double(bnd_X.(x{1}) ~= 0 | bnd_X.(x{2}) ~= 0);
            h3.AlphaDataMapping = 'none';
            h3.FaceColor = 'flat';
            h3.EdgeColor = 'none';
            h3.AlphaData = flow_boundary;
            h3.FaceAlpha = 'flat';
        end
    end
    
    % Modify viewing
    ax.View = [0 90];
    ax.Box = 'on';
    ax.NextPlot = 'replacechildren';
    colormap(ax, 'jet');
    axis(ax, 'equal')
    axis(ax, 'tight');
    axis(ax, 'manual');
else
    % If plot already present update values
    h_mag.CData = data.field;
end
drawnow;

% Return handles
if nargout == 1
    h_magnitude = h_mag;
end
if nargout == 2
    h_magnitude = h_mag;
    cax = ax;
end