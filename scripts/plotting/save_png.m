files = dir();
for i = 1:length(files)
    if files(i).isdir
        continue;
    end
    h = hgload(files(i).name);
    name = files(i).name;
    pattern = '.fig';
    replace = '';
    name = regexprep(name, pattern, replace);
    saveas(h, strrep(name, ' ', '_'), 'png')
end