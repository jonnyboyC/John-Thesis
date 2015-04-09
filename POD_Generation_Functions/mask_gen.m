function [bnd_x, bnd_y] = mask_gen(data, bnd_idx, bnd_x, bnd_y)
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','on','Position',[500, 500, 900, 500]);

% Construct the components.
hupdate = uicontrol('Style','pushbutton',...
             'String','Update Plot','Position',[500,50,100,50],...
             'Callback',@hupdate_Callback);
hsave   = uicontrol('Style','pushbutton',...
             'String','Save / Exit','Position',[760,50,100,50],...
             'Callback',@hsave_Callback); 
hclear   = uicontrol('Style','pushbutton', ...
             'String','Clear','Position',[630,50,100,50],...
             'Callback',@hclear_Callback);
         
hcheck1     = uicontrol('Style', 'checkbox', ...
                'String', 'Left', ...
                'Value', 1, ...
                'Position', [500, 360, 60, 20]);
hcheck2     = uicontrol('Style', 'checkbox', ...
                'String', 'Right', ...
                'Value', 1, ...
                'Position', [580, 360, 60, 20]);
hcheck3     = uicontrol('Style', 'checkbox', ...
                'String', 'Bottom', ...
                'Value', 1, ...
                'Position', [660, 360, 60, 20]);
hcheck4     = uicontrol('Style', 'checkbox', ...
                'String', 'Top', ...
                'Value', 1, ... 
                'Position', [740, 360, 60, 20]);
    
htext1   = uicontrol('Style', 'text', ...
                    'String', 'Exclusion Zones', ...
                    'FontSize', 16, ...
                    'Position', [600, 400, 160, 50]);
htext2   = uicontrol('Style', 'text', ...
                    'String', ['Boundaries Image ' num2str(size(data.x,1)) ' x ' num2str(size(data.x,2))], ...
                    'FontSize', 16, ...
                    'Position', [100, 400, 300, 50]);
         
columnname = {'Left', 'Right', 'Bottom', 'Top'};
table_data = cell(8, 4);

htable = uitable(f, 'Data', table_data, ...
                    'ColumnName', columnname, ...
                    'ColumnEditable', true(1,4), ...
                    'RowName', [], ...
                    'Position', [500, 150, 310, 200]);
haxes = axes('Units','pixels','Position',[50,50,400,400]);

% Initialize the GUI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hsave.Units = 'normalized';
hupdate.Units = 'normalized';
hclear.Units = 'normalized';
htable.Units = 'normalized';
hcheck1.Units = 'normalized';
hcheck2.Units = 'normalized';
hcheck3.Units = 'normalized';
hcheck4.Units = 'normalized';
htext1.Units = 'normalized';
htext2.Units = 'normalized';

bnd_x_temp = bnd_x;
bnd_y_temp = bnd_y;
% Generate the data to plot.
[~, haxes] = Plottec2(data, haxes);

% Assign the GUI a name to appear in the window title.
f.Name = 'Open Flow Boundary Refinement';

% Move the GUI to the center of the screen.
movegui(f,'center')

% Make the GUI visible.
uiwait(f);
%  Pop-up menu callback. Read the pop-up menu Value property to
%  determine which item is currently displayed and make it the
%  current data. This callback automatically has access to 
%  current_data because this function is nested at a lower level.

  % Push button callbacks. Each callback plots current_data in the
  % specified plot type.

    function hsave_Callback(~, ~) 
    % Display mesh plot of the currently selected data.
        bnd_x = bnd_x_temp;
        bnd_y = bnd_y_temp;
        uiresume(f);
        close all
    end

    function hupdate_Callback(~, ~) 
    % Display contour plot of the currently selected data.
        bounds = htable.Data;

        bnd_x_temp = bnd_x;
        bnd_y_temp = bnd_y;
        if ~hcheck1.Value
            bnd_y_temp(bnd_y_temp > 0) = 0;
        end
        if ~hcheck2.Value
            bnd_y_temp(bnd_y_temp < 0) = 0;
        end
        if ~hcheck3.Value
            bnd_x_temp(bnd_x_temp > 0) = 0;
        end
        if ~hcheck4.Value
            bnd_x_temp(bnd_x_temp < 0) = 0;
        end
        if ~all(cellfun('isempty', bounds))
            filled = logical(all(~cellfun('isempty', bounds)')'*ones(1, 4));
            filters = cellfun(@str2double, bounds(filled));
            filters = reshape(filters, [], 4);
            for i = 1:size(filters,1);
                if filters(i,1) <= 0
                    filters(i,1) = 1;
                end
                if filters(i,2) > size(data.x,1)
                    filters(i,2) = size(data.x,1);
                end
                if filters(i,1) >= filters(i,2)
                    filters(i,2) = filters(i,1) + 1;
                end
                if filters(i,3) <= 0
                    filters(i,3) = 1;
                end
                if filters(i,4) > size(data.x,2)
                    filters(i,4) = size(data.x,2);
                end
                if filters(i,3) >= filters(i,4)
                    filters(i,4) = filters(i,3) + 1;
                end
                bnd_x_temp(filters(i,1):filters(i,2), filters(i,3):filters(i,4)) = 0;
                bnd_y_temp(filters(i,1):filters(i,2), filters(i,3):filters(i,4)) = 0;
            end
            htable.Data(1:size(filters,1),:) = ...
                cellfun(@num2str, num2cell((filters)), 'UniformOutput', false);
        end
        data.bnd_x = bnd_x_temp;
        data.bnd_y = bnd_y_temp;
        [~, haxes] = Plottec2(data, haxes);
    end

    function hclear_Callback(~, ~)
        htable.Data = cell(8, 4);
        data.bnd_x = bnd_x;
        data.bnd_y = bnd_y;
        [~, haxes] = Plottec2(data, haxes);
    end

end