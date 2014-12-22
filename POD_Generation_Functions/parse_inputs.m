function problem = parse_inputs(fields, setdefaults, varargin)
% Parse inputs into correct stuture, using passed cell array of FIELDS to
% specify field names and funciton handle SETDEFAULTS to set correct defaults 


% if no input is provided or a variable create a blank structure
if nargin == 2 || ~isstruct(varargin{1})
    fprintf('Non-Struct provided, setting POD_Gen to defaults\n\n');
    problem = struct(); 
else
    problem = varargin{1};
end

% Find fields that have already been set
curr_fields = isfield(problem, fields);

% create fields that are not present
for i = 1:size(curr_fields, 2)
    if ~curr_fields(i)
       problem.(fields{i}) = [];
    end
end

% Find any extra fields that have been provided that don't belong
field_names = fieldnames(problem);
rm_fields = field_names(~ismember(field_names, fields));

% Remove fields and notify user
if ~isempty(rm_fields)
    for i = 1:size(rm_fields,1);
        fprintf('Removing field "%s"\n', rm_fields{i});
    end
    problem = rmfield(problem, rm_fields);
end

problem = setdefaults(problem);
end