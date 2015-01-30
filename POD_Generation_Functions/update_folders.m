function update_folders(direct)
%% Update data folders to follow correct formatting

% Check if the data folders are present, create them if they are not.
% These include a folder for raw vc7/similar files, a processed vc7
% data, folder post pod calculation, and a figures folder
if ~exist([direct '\Processed Data'], 'dir')
    mkdir(direct, '\Processed Data');
end
if ~exist([direct '\Other'], 'dir')
    mkdir(direct, '\Other');
end
if ~exist([direct '\Raw Data'], 'dir')
    mkdir(direct, '\Raw Data');
end
if ~exist([direct '\POD Data'], 'dir')
    mkdir(direct, '\POD Data');
end
if ~exist([direct '\Galerkin Coeff'], 'dir')
    mkdir(direct, '\Galerkin Coeff');
end
if ~exist([direct '\Mod Galerkin Coeff'], 'dir')
    mkdir(direct, '\Mod Galerkin Coeff')
end
if ~exist([direct '\Viscous Coeff'], 'dir')
    mkdir(direct, '\Viscous Coeff');
end
if ~exist([direct '\Figures'], 'dir')
   mkdir(direct, '\Figures');
   if ~exist([direct '\Figures\POD'], 'dir')
       mkdir(direct, '\Figures\POD');
   end
   if ~exist([direct '\Figures\Galerkin'], 'dir')
       mkdir(direct, '\Figures\Galerkin');
   end  
   if ~exist([direct '\Figures\Movies'], 'dir')
       mkdir(direct, '\Figures\Movies');
   end 
end
% Move any vc7 or mat files into the Raw data folder
[~, ~, ~] = movefile([direct '\*.vc7'], [direct '\Raw Data']); 
[~, ~, ~] = movefile([direct '\*.im7'], [direct '\Raw Data']); 
[~, ~, ~] = movefile([direct '\*.mat'], [direct '\Raw Data']); 
[~, ~, ~] = movefile([direct '\*.dat'], [direct '\Raw Data']); 
end