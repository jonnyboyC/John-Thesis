close all
files = dir();
for i = 1:length(files)
    if files(i).isdir
        continue;
    end
    h = hgload(files(i).name);
    ax = gca;

    % FORMAT
    cluster_num_format(ax)
    % FORMAT
    
    name = files(i).name;
    pattern = '.fig';
    replace = '';
    name = regexprep(name, pattern, replace);
    saveas(h, strrep(name, ' ', '_'), 'png')
    saveas(h, name, 'fig');
end
