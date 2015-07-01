function [bnd_X, bnd_idx] = mask_gen(data, U, X_direct, streamlines)
% MASK_GEN generate a mask for the raw PIV data
% [bnd_x, bnd_y, bnd_idx] = MASK_GEN(data, bnd_x, bnd_y, u, v) open a 
% GUI to exclude sections of the flow images that may be causing problems
% as well as correctly identifying where open flow boundaries are occuring

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','on','Position',[500, 500, 900, 500]);

% Construct control buttons
hupdate = uicontrol('Style','pushbutton',...
             'String','Update Plot','Position',[500,50,100,50],...
             'Callback',@hupdate_Callback);
hsave   = uicontrol('Style','pushbutton',...
             'String','Save / Exit','Position',[760,50,100,50],...
             'Callback',@hsave_Callback); 
hclear   = uicontrol('Style','pushbutton', ...
             'String','Clear','Position',[630,50,100,50],...
             'Callback',@hclear_Callback);
         
% Construct checkboxes
hcheck1     = uicontrol('Style', 'checkbox', ...
                'String', 'Left', ...
                'Value', 1, ...
                'Position', [500, 420, 60, 20]);
hcheck2     = uicontrol('Style', 'checkbox', ...
                'String', 'Right', ...
                'Value', 1, ...
                'Position', [580, 420, 60, 20]);
hcheck3     = uicontrol('Style', 'checkbox', ...
                'String', 'Bottom', ...
                'Value', 1, ...
                'Position', [660, 420, 60, 20]);
hcheck4     = uicontrol('Style', 'checkbox', ...
                'String', 'Top', ...
                'Value', 1, ... 
                'Position', [740, 420, 60, 20]);
   
x = flow_comps_ns(data.X);
            
% Construct textboxes
htext1   = uicontrol('Style', 'text', ...
                    'String', 'Flow Boundary Exclusion Zones', ...
                    'FontSize', 16, ...
                    'Position', [570, 440, 160, 50]);
htext2   = uicontrol('Style', 'text', ...
                    'String', 'Flow Exclusion Zones', ...
                    'FontSize', 16, ...
                    'Position', [570, 250, 160, 50]);
htext3   = uicontrol('Style', 'text', ...
                    'String', ['Boundaries Image ' num2str(size(data.X.(x{1}),1)) ' x ' num2str(size(data.X.(x{1}),2))], ...
                    'FontSize', 16, ...
                    'Position', [80, 440, 300, 50]);
         
% Produce Tables
columnname = {'Left', 'Right', 'Bottom', 'Top'};
table_data = cell(4, 4);

htable1 = uitable(f, 'Data', table_data, ...
                    'ColumnName', columnname, ...
                    'ColumnEditable', true(1,4), ...
                    'RowName', [], ...
                    'Position', [500, 310, 310, 100]);
                
htable2 = uitable(f, 'Data', table_data, ...
                    'ColumnName', columnname, ...
                    'ColumnEditable', true(1,4), ...
                    'RowName', [], ...
                    'Position', [500, 130, 310, 100]);
axes('Units','pixels','Position',[50,50,400,400]);

% Initialize the GUI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hsave.Units = 'normalized';
hupdate.Units = 'normalized';
hclear.Units = 'normalized';
htable1.Units = 'normalized';
htable2.Units = 'normalized';
hcheck1.Units = 'normalized';
hcheck2.Units = 'normalized';
hcheck3.Units = 'normalized';
hcheck4.Units = 'normalized';
htext1.Units = 'normalized';
htext2.Units = 'normalized';
htext3.Units = 'normalized';

% Create a temporary copy of the intially calculated values
bnd_X_temp = data.bnd_X;
bnd_idx_temp = data.bnd_idx;

% Generate a plot of the mean flow
plot_vector_field(data, streamlines);
colorbar;

[u, x] = flow_comps(U, data.X);
dims = flow_dims(U);
images = size(U.(u{1}), dims);
dimensions = size(data.X.(x{1}));

% Assign the GUI a name to appear in the window title.
f.Name = 'Open Flow Boundary Refinement';

% Move the GUI to the center of the screen.
movegui(f,'center')

% Halt execution of program to allow for user input
uiwait(f);

% Nested callback functions to control the plot
    function hsave_Callback(~, ~)
        % Save values and resume execution
        
        % Set the temp values as the new values
        bnd_X = bnd_X_temp;
        bnd_idx = bnd_idx_temp;
        
        % Resume execution and close plot
        uiresume(f);
        close(f)
    end

    function hupdate_Callback(~, ~) 
        % Update plot with values from table
        
        % Pull data from tables
        bnd_bounds = htable1.Data;
        flow_bounds = htable2.Data;
        
        % Intially set all to in flow
        bnd_idx_temp = ones(dimensions);
                    
        idx = struct_index({[1 1]}, dims(end), U);
        comps = flow_ncomps(U);
        
        % Make all images with a pixel out of flow as part of boundary
        for i = 1:images
            temp = zeros(dimensions);
            for j = 1:comps
                idx{end} = i;
                temp = temp + U.(u{j})(idx{:});
            end
            bnd_idx_temp = bnd_idx_temp + double(temp == 0);
        end
        
        % Set all points with at least one image not captured as in boundary
        bnd_idx_temp(bnd_idx_temp > images/100) = -1;
        bnd_idx_temp(bnd_idx_temp > 1) = 1;
        

        % Check that a full row has been filled in 
        if any(~all(cellfun('isempty', flow_bounds)'))
            
            % Pull filters from table
            filled = logical(all(~cellfun('isempty', flow_bounds)')'*ones(1, 4));
            filters = cellfun(@str2double, flow_bounds(filled));
            filters = reshape(filters, [], 4);
            
            % Apply each filter to the data
            for i = 1:size(filters,1);
                if filters(i,1) <= 0
                    filters(i,1) = 1;
                end
                if filters(i,2) > dimensions(1)
                    filters(i,2) = dimensions(1);
                end
                if filters(i,1) >= filters(i,2)
                    filters(i,2) = filters(i,1) + 1;
                end
                if filters(i,3) <= 0
                    filters(i,3) = 1;
                end
                if filters(i,4) > dimensions(2)
                    filters(i,4) = dimensions(2);
                end
                if filters(i,3) >= filters(i,4)
                    filters(i,4) = filters(i,3) + 1;
                end
                % Set filtered areas to -1 for in boundary
                bnd_idx_temp(filters(i,1):filters(i,2), filters(i,3):filters(i,4)) = -1;
            end
            htable2.Data(1:size(filters,1),:) = ...
                cellfun(@num2str, num2cell((filters)), 'UniformOutput', false);
        end
        
        if ~open_flow
            % Use built in edge detection to boundary points change points to 0
            bnd_idx_temp(edge(bnd_idx_temp, 'sobel')) = 0;
            
            % manually change edge points on image edge
            bnd_idx_temp = manual_edge(bnd_idx_temp);
        end
        
        % Determine potential flow boundaries
        bnd_X_temp = edge_boundaries(bnd_idx_temp, data.X, X_direct);
        
        % Apply checkboxes
        if ~hcheck1.Value
            bnd_X_temp.(x{1})(bnd_X_temp.(x{1}) > 0) = 0;
        end
        if ~hcheck2.Value
            bnd_X_temp.(x{1})(bnd_X_temp.(x{1}) < 0) = 0;
        end
        if ~hcheck3.Value
            bnd_X_temp.(x{2})(bnd_X_temp.(x{2}) > 0) = 0;
        end
        if ~hcheck4.Value
            bnd_X_temp.(x{2})(bnd_X_temp.(x{2}) > 0) = 0;
        end
        
        % Check that a full row has been filled in 
        if any(~all(cellfun('isempty', bnd_bounds)'))
            
            % Pull filters from table
            filled = logical(all(~cellfun('isempty', bnd_bounds)')'*ones(1, 4));
            filters = cellfun(@str2double, bnd_bounds(filled));
            filters = reshape(filters, [], 4);
            
            % Apply each filter to the data
            for i = 1:size(filters,1);
                if filters(i,1) <= 0
                    filters(i,1) = 1;
                end
                if filters(i,2) > dimensions(1)
                    filters(i,2) = dimensions(1);
                end
                if filters(i,1) >= filters(i,2)
                    filters(i,2) = filters(i,1) + 1;
                end
                if filters(i,3) <= 0
                    filters(i,3) = 1;
                end
                if filters(i,4) > dimensions(2)
                    filters(i,4) = dimensions(2);
                end
                if filters(i,3) >= filters(i,4)
                    filters(i,4) = filters(i,3) + 1;
                end
                % Set filtered areas to 0 to indicate it is not an open
                % flow boundary
                bnd_X_temp.(x{1})(filters(i,1):filters(i,2), filters(i,3):filters(i,4)) = 0;
                bnd_X_temp.(x{1})(filters(i,1):filters(i,2), filters(i,3):filters(i,4)) = 0;
            end
            htable1.Data(1:size(filters,1),:) = ...
                cellfun(@num2str, num2cell((filters)), 'UniformOutput', false);
        end
        data.bnd_X = bnd_X_temp;
        data.bnd_idx = bnd_idx_temp;
        
        % Re-plot
        plot_vector_field(data, streamlines);
        colorbar;
    end

    function hclear_Callback(~, ~)
        % Clear values from table 
        htable1.Data = cell(4, 4);
        htable2.Data = cell(4, 4);
        data.bnd_X = data.bnd_X;
        data.bnd_idx = data.bnd_idx;
        
        % Re-plot
        plot_vector_field(data, streamlines);
        colorbar;
    end

end