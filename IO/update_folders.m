function update_folders(direct)
%% Update data folders to follow correct formatting

% Check if the data folders are present, create them if they are not.
% These include a folder for raw vc7/similar files, a processed vc7
% data, folder post pod calculation, and a figures folder
if ~exist([direct filesep 'Processed Data'], 'dir')
    mkdir([direct filesep 'Processed Data']);
end
if ~exist([direct filesep 'Processed Data' filesep 'Mask'], 'dir')
   mkdir([direct filesep 'Processed Data' filesep 'Mask']);
end
if ~exist([direct filesep 'Other'], 'dir')
    mkdir([direct filesep 'Other']);
end
if ~exist([direct filesep 'Raw Data'], 'dir')
    mkdir([direct filesep 'Raw Data']);
end
if ~exist([direct filesep 'POD Data'], 'dir')
    mkdir([direct filesep 'POD Data']);
end
if ~exist([direct filesep 'Galerkin Coeff'], 'dir')
    mkdir([direct filesep 'Galerkin Coeff']);
end
if ~exist([direct filesep 'Mod Galerkin Coeff'], 'dir')
    mkdir([direct filesep 'Mod Galerkin Coeff'])
end
if ~exist([direct filesep 'Viscous Coeff'], 'dir')
    mkdir([direct filesep 'Viscous Coeff']);
end
if ~exist([direct filesep 'Figures'], 'dir')
   mkdir([direct filesep 'Figures']);
end
if ~exist([direct filesep 'Figures' filesep 'POD'], 'dir')
   mkdir([direct filesep 'Figures' filesep 'POD']);
end
if ~exist([direct filesep 'Figures' filesep 'POD' filesep 'Modes'], 'dir')
   mkdir([direct filesep 'Figures' filesep 'POD' filesep 'Modes']);
end
if ~exist([direct filesep 'Figures' filesep 'POD' filesep 'Clusters'], 'dir')
   mkdir([direct filesep 'Figures' filesep 'POD' filesep 'Clusters']);
end
if ~exist([direct filesep 'Figures' filesep 'Galerkin'], 'dir')
   mkdir([direct filesep 'Figures' filesep 'Galerkin']);
end  
if ~exist([direct filesep 'Figures' filesep 'Movies'], 'dir')
   mkdir([direct filesep 'Figures' filesep 'Movies']);
end 
% Move any vc7 or mat files into the Raw data folder
[~, ~, ~] = movefile([direct filesep '*.vc7'], [direct filesep 'Raw Data']); 
[~, ~, ~] = movefile([direct filesep '*.VC7'], [direct filesep 'Raw Data']); 
[~, ~, ~] = movefile([direct filesep '*.im7'], [direct filesep 'Raw Data']); 
[~, ~, ~] = movefile([direct filesep '*.IM7'], [direct filesep 'Raw Data']); 
[~, ~, ~] = movefile([direct filesep '*.mat'], [direct filesep 'Raw Data']); 
[~, ~, ~] = movefile([direct filesep '*.dat'], [direct filesep 'Raw Data']); 

% TODO want to select all other file types that arent' folders
[~, ~, ~] = movefile([direct filesep '*.xls'], [direct filesep 'Other']);
end