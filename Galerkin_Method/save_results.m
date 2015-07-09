function save_results(num_modes, direct, folder, custom, varargin)
% check if folder exist create if empty
if custom
    direct_ext = [direct filesep folder filesep 'modes_' num_modes '_custom'];
else
    direct_ext = [direct filesep folder filesep 'modes_' num2str(num_modes)];
end

if ~exist(direct_ext, 'dir')
    mkdir(direct_ext);
end

% file name for mat file
filename = [direct_ext filesep 'Coefficients_run_' num2str(varargin{1}.run_num) '.mat'];

% allow direct access to file
file = matfile(filename, 'Writable', true);

% save all 
for i = 1:length(varargin)
    file.(varargin{i}.name) = varargin{i};
end
end