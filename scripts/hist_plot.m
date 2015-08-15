function hist_plot(big_names, big_R, varargin)

direct = 'C:\Users\John-Desktop\Desktop\histograms\Cavity';

sub_list = zeros(length(big_R));
switch nargin
    case 3
        tag1 = varargin{1};
        
        for i = 1:length(big_R)
            if ~isempty(strfind(big_names{i}, tag1)) 
                sub_list(i) = big_R(i);
            end
        end
    case 4
        tag1 = varargin{1};
        tag2 = varargin{2};
        
        for i = 1:length(big_R)
            if (~isempty(strfind(big_names{i}, tag1)) && ~isempty(strfind(big_names{i}, tag2)))
                sub_list(i) = big_R(i);
            end
        end
    case 5
        tag1 = varargin{1};
        tag2 = varargin{2};
        tag3 = varargin{3};
        
        for i = 1:length(big_R)
            if (~isempty(strfind(big_names{i}, tag1)) ...
                    && ~isempty(strfind(big_names{i}, tag2)) ...
                    && ~isempty(strfind(big_names{i}, tag3)))
                sub_list(i) = big_R(i);
            end
        end
end


        
% Yeah none of these happen
sub_list(sub_list == 0) = [];

vert_lines = mean(sub_list);

h = figure;
hist(sub_list, 20);
hold on
ax = gca;
plot([vert_lines, vert_lines;], [ax.YLim(2) ax.YLim(1)], '-.r')
legend('Correlation bins', ['Mean Correlation ' num2str(vert_lines)], 'Location', 'northwest');

switch nargin
    case 3
        saveas(h, [direct filesep tag1], 'fig');
    case 4
        saveas(h, [direct filesep tag1 '_v_' tag2], 'fig');
    case 5
        saveas(h, [direct filesep tag1 '_v_' tag2 '_v_' tag3], 'fig');
end

ax.XLabel = xlabel('Correlation Value');
ax.YLabel = ylabel('Correlations');
ax.XLim = [-1 1];


end