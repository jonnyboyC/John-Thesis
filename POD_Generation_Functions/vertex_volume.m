function volume = vertex_volume(X, bnd_idx)

% Get dimensions information
dims = flow_dims(X);
x = flow_comps_ip(X);
dimensions = size(bnd_idx);

% Set all points for open flow on edge to on boundary 
exterior = false(size(bnd_idx));
for i = 1:dims
    idx_1 = flow_index([1, 1], i, exterior);
    idx_end = flow_index([0, 0], i, exterior);
    
    exterior(idx_1{:}) = true;
    exterior(idx_end{:}) = true;
end
bnd_idx(exterior & bnd_idx == 1) = 0;

% Wrap exterior with in boundary 
bnd_temp = bnd_idx;
bnd_idx = -ones(dimensions+2);
idx = flow_index(repmat({[2 -1]}, 1, dims), 1:dims, bnd_idx);
bnd_idx(idx{:}) = bnd_temp;

% Get choices for shift
median_shifts{i} = cell(1,dims);
for i = 1:dims
    median_shifts{i} = get_shifts(i);
end

% All interstital points on the mesh for an arbitrary mesh
X_median = [];
for i = 1:dims
    idx = flow_index(repmat({[1 -1]}, 1, dims), 1:dims, X.(x{i}));
    X_temp = [];
    for j = 1:length(median_shifts)  
        for k = 1:length(median_shifts{j})
            
            % determine all combinations of dimension choices
            averaged_dims = nchoosek(1:dims, size(median_shifts{j},2));
            
            for l = 1:size(averaged_dims,1)
                
                comp = '';
                idx_temp = idx;
                cnt = 1;
                
                % Alter index to create proper summing
                for ll = averaged_dims(l,:)
                    idx_temp{ll} = idx_temp{ll} + median_shifts{j}(k, cnt);
                    comp = [comp x{ll}];
                    cnt = cnt + 1;
                end
                
                % Add to previous results
                if isfield(X_temp, comp);
                    X_temp.(comp) = X_temp.(comp) ...
                        + X.(x{i})(idx_temp{:})/length(median_shifts{j});
                else
                    X_temp.(comp) = X.(x{i})(idx_temp{:})/length(median_shifts{j});
                end
            end
        end
    end
    X_median.(x{i}) = X_temp;
end

% Preallocate values
dimensions = size(bnd_idx);
volume = zeros(dimensions);
idx = cell(1, dims);
idx_temp = cell(1, dims);
mask_idx = flow_index(repmat({[1 2]}, 1, dims), 1:dims, bnd_idx);
shifts = get_shifts(dims);
xy = flow_comps(X_median.(x{1}));


% Determine volume
for i = 1:numel(bnd_idx)
    
    % If in boundary set to zero and continue, includes exterior buffer
    if bnd_idx(i) == -1;
        volume(i) = 0;
        continue;
    end
    
    % Get n-dimensional index
    [idx{:}] = ind2sub(dimensions, i);
    
    % If open flow get volume as box surround vertex, else calculate
    % partial box for on boundary points
    if bnd_idx(i) == 1;
        mask_temp = mask_idx;
        vertices = zeros(pow2(dims),dims);
        for j = 1:dims
            mask_temp{j} = mask_temp{j} + idx{j} - 3;
        end
        for j = 1:dims
            vertices(:,j) = reshape(X_median.(x{j}).(sIdx(x, 1:dims))(mask_temp{:}),[],1);
        end
        [~, volume(i)] = convhull(vertices);
    else
        for j = 1:size(shifts,1)
            mask_temp = mask_idx;
            for k = 1:dims
                mask_temp{k} = mask_temp{k} + idx{k} + shifts(j,k) - 2;
            end
            
            % If we find a corner with all on boundary or in flow include corner
            if all(bnd_idx(mask_temp{:}) >= 0)
                vertices = zeros(pow2(dims),dims);
                for k = 1:dims
                    vertices(1,k) = X.(x{k})(i);
                end
                for k = 1:dims
                    for l = 1:dims
                        idx_temp = num2cell([idx{:}] - 1);
                        idx_temp{k} = idx_temp{k} + shifts(j,k) - 1;
                        vertices(l+1,k) = X_median.(x{k}).(xy{l})(idx_temp{:});
                    end
                end
                [~, vol_temp] = convhull(vertices);
                volume(i) = volume(i) + vol_temp;
            end
        end
    end
end

idx = flow_index(repmat({[2 -1]}, 1, dims), 1:dims, volume);
volume = volume(idx{:});

end

function idx = sIdx(x, range)
% Return string to the compound structure name

idx = [];
for i = range
   idx = [idx x{i}]; 
end
end