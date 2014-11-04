function run_all

start_direct = 'D:\shear layer';
direct = uigetdir(start_direct, 'Choose Source Image Directory');
folders = dir([direct '\*']);

for i = 3:length(folders)
   	if folders(i).isdir
       main_code_chabot(20, true, true, false, 'jpg', [direct '\' folders(i).name])
    end
end