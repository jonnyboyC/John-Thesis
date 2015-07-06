function transform = inverse_jacobian(XdXi, dimensions)

% Get flow information
xi = flow_comps(XdXi);
xi_comps = flow_ncomps(XdXi);
elems = prod(dimensions);

% Preallocat
jacobian = zeros(xi_comps, xi_comps, elems);
jacobian_inv = zeros(xi_comps, xi_comps, elems);

% Create jacobain matrix for each derivative component
for i = 1:xi_comps
    for j = 1:xi_comps
       jacobian(i,j,:) = XdXi.(xi{i}).(xi{j})(:);
    end
end

% Take inverse jacobian using matlab's inv for more than 2 dimensions
if length(dimensions) == 2
    for i = 1:elems
        jacobian_inv(:,:,i) = 1/det(jacobian(:,:,i))*[jacobian(2,2), -jacobian(1,2);
                                                      -jacobian(2,1), jacobian(1,1)];
    end
else
    for i = 1:elems
        jacobian_inv(:,:,i) = inv(jacobian(:,:,i));
    end
end

% Convert back into structured form
for i = 1:xi_comps
    for j = 1:xi_comps
        transform.(xi{i}).(xi{j}) = reshape(jacobian_inv(i,j,:), dimensions);
        transform.(xi{i}).(xi{j})(isinf(transform.(xi{i}).(xi{j}))) = 0;
    end
end
end