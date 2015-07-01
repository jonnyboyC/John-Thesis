function shifts = get_shifts(dims)
% temporary super crummy way to do this

choices = [zeros(1, dims), ones(1, dims)];
temp = unique(nchoosek(choices, dims), 'rows');

shifts = [];
for i = 1:size(temp,1)
     shifts = [shifts; unique(perms(temp(i,:)), 'rows')];
end
end