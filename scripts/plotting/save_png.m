files = dir();
for i = 1:length(files)
    if files(i).isdir
        continue;
    end
    h = hgload(files(i).name);
    % REMOVE
    ax = gca;
    ax.Children(2).CData = ax.Children(2).CData*0.3048
    % REMOVE
    name = files(i).name;
    pattern = '.fig';
    replace = '';
    name = regexprep(name, pattern, replace);
    saveas(h, strrep(name, ' ', '_'), 'png')
end