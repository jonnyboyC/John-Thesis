function ao = til_init(modal_amp, init, RD_nm, X) 
% basis transfrom for inital conditions
ao = zeros(1,RD_nm);
for i = 1:RD_nm
    ao(1, i) = sum(X(:,i)'.*modal_amp(init,:));
end