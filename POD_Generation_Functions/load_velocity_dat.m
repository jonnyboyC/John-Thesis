function [u, v, w] = load_velocity_dat(direct, num_images num_zones)
% TODO make this not the worst code ever, why are we using random numbers
% to read a file??
 %#ok<*AGROW>

u = [];
v = [];
w = [];

u_zones = cells(num_zones, 1);
v_zones = cells(num_zones, 1);
w_zones = cells(num_zones, 1);

velocity_file = fopen([direct filesep 'Raw Data' filesep 'velocity.dat']);
i = 'integer*4';
f = 'single';

% May need to change
skips = 1;
image_count = 1;

for image = 1:num_images
    % I believe this has to do with selectively not sampling 
    if abs(round((image-1)/skips)-((image-1)/skips)) > 0
        for zone = 1:num_zones
            % Get image dimensions and determine the number of
            % positions to read in the file
            dimensions = [fread(velocity_file,1,i), fread(velocity_file,1,i)];
            read_length = prod(dimensions);
            
            % Appears we are reading maybe randomly throwing away some
            % data????
            fread(velocity_file,read_length,f);
            fread(velocity_file,read_length,f);
            fread(velocity_file,read_length,f);
            fread(velocity_file,read_length,f);
        end
    else
        for zone = 1:num_zones
            % Get image dimensions, probably the same in each will most
            % likely change
            dimensions = [fread(velocity_file,1,i), fread(velocity_file,1,i)];
            
            % Get number of intergers to read based on dimensions
            read_length = prod(dimensions);
            
            % Read unneeded radial velocity component
            fread(velocity_file,read_length,f);
            
            % Original file performed some scaling, not sure if this is
            % a combined scaling and unit conversion
            u_temp = reshape(fread(velocity_file,read_length,f), dimensions);
            v_temp = reshape(fread(velocity_file,read_length,f), dimensions);
            w_temp = reshape(fread(velocity_file,read_length,f), dimensions);
            
            u_zones{zone}(:,:,image_count) = u_temp;
            v_zones{zone}(:,:,image_count) = v_temp;
            w_zones{zone}(:,:,image_count) = w_temp;
        end
        image_count = image_count + 1;
    end
end

fclose(velocity_file);

% Stack data hopefully the correct way
for zone = 1:num_zones
    u = [u; u_zones{zone}];
    v = [v; v_zones{zone}];  
    w = [w; w_zones{zone}];
end
end