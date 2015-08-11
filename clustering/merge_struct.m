function [target, name] = merge_struct(var, stop_on)
if ~isstruct(var)
    target = [];
    name = {};
    return;
end

fields = fieldnames(var);
if any(ismember(fields, stop_on))
    match = ismember(fields, stop_on);
    if var.completed 
        target = var.(fields{match});
        name = {};
    else
        target = [];
        name = {};
    end
    return;
end

models = {'GM', 'GM1', 'GM2', 'GM3'};
target = [];
name = {};
for i = 1:length(fields)
    [target_nest, name_nest] = merge_struct(var.(fields{i}), stop_on);
    if any(ismember(fields, models));
        for j = 1:length(name_nest)
            name_nest{j} = [fields{i} '_' name_nest{j}];
        end
        name = [name, name_nest];
    elseif isscalar(target_nest)
        name = [name, fields(i)];
    elseif ~isempty(name_nest)
        name = [name, name_nest];
    end
    target = [target, target_nest];
end