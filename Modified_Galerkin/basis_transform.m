function [pod_U_til, modal_amp_til] = ...
                basis_transform(pod_Ut, modal_amp, num_modes, X) 
% Transform pod modes and modal amplitude into new basis
comps = flow_ncomps(pod_Ut);
u = flow_comps(pod_Ut);
for i = 1:comps
    pod_U_til.(u{i}) = zeros(size(pod_Ut.(u{i}),1), num_modes);
end

modal_amp_til = zeros(size(modal_amp,1), num_modes);

for i = 1:num_modes
    for j = 2:size(pod_Ut.(u{1}),1)
        for k = 1:comps
            pod_U_til.(u{k})(j,i) = sum(X(:,i)'.*pod_Ut.(u{k})(j,:));
        end
    end
    for j = 2:size(modal_amp,1)
        modal_amp_til(j,i) = sum(X(:,i)'.*modal_amp(j,:));
    end
end

end