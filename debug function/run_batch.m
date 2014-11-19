function run_batch

% main_code_chabot(2000, false, true, false, 'jpg', 'D:\shear layer\PIVData\Forced1');
% main_code_chabot(2000, false, true, false, 'jpg', 'D:\shear layer\PIVData\Baseline1');

gcp;
% parfor i = 1:9
%     a = 4:2:20;
%     Galerkin_Proj(a(i), {'amp', 'fft'}, true, 0:0.0001:50, 1, 'D:\shear layer\PIVData\Baseline1');
% end

parfor i = 1:10
    a = 2:2:20;
    Galerkin_Proj(a(i), 'amp', true, 0:0.0001:50, 1, 'D:\shear layer\PIVData\Forced1')
end
% end

% Mod_POD(4, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Baseline1', 'Coeff_m8i1.mat');
% Mod_POD(6, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Baseline1', 'Coeff_m12i1.mat');
% Mod_POD(8, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Baseline1', 'Coeff_m16i1.mat');
% Mod_POD(10, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Baseline1', 'Coeff_m20i1.mat');
% Mod_POD(12, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Baseline1', 'Coeff_m24i1.mat');

Mod_POD(4, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Forced1', 'Coeff_m8i1.mat');
Mod_POD(6, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Forced1', 'Coeff_m12i1.mat');
Mod_POD(8, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Forced1', 'Coeff_m16i1.mat');
Mod_POD(10, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Forced1', 'Coeff_m20i1.mat');
Mod_POD(12, {'amp', 'fft'}, true, 1, 100, 'D:\shear layer\PIVData\Forced1', 'Coeff_m24i1.mat');