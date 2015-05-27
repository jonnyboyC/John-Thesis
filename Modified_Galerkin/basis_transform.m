function [pod_u_til, pod_v_til, modal_amp_til] = ...
                basis_transform(pod_ut, pod_vt, modal_amp, RD_nm, X) 
% Transform pod modes and modal amplitude into new basis
pod_u_til = zeros(size(pod_ut,1), RD_nm);
pod_v_til = zeros(size(pod_vt,1), RD_nm);
modal_amp_til = zeros(size(modal_amp,1), RD_nm);
for i = 1:RD_nm
    for j = 2:size(pod_ut,1)
        pod_u_til(j,i) = sum(X(:,i)'.*pod_ut(j,:));
        pod_v_til(j,i) = sum(X(:,i)'.*pod_vt(j,:));
    end
    for j = 2:size(modal_amp,1)
        modal_amp_til(j,i) = sum(X(:,i)'.*modal_amp(j,:));
    end
end

end