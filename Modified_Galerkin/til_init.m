function ao = til_init(modal_amp, init, num_modes, X) 
% TIL_INIT apply basis trasnform to initial conditions of intergration
%
%   ao = TIL_INIT(modal_amp, init, 

ao = zeros(1,num_modes);
for i = 1:num_modes
    ao(1, i) = sum(X(:,i)'.*modal_amp(init,:));
end