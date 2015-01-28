function u_scale = u_scale_gen_airfoil(~, direct)
airfoil_path = [direct, 'other\A2_CpPHF_NS_20100728b.xls';
R = 287;
water2Pa = 249.17;
k = .952;  % Tunnel constant
patm = 98600; % Pa
data = xlsread(airfoil_path, 'DataSheet'); % read data from excel sheet, "act" not used here
dp = water2Pa*data(:,4)/k;  % Dynamic Pressure (Pa)
sp = water2Pa*data(:,5);    % Static Pressure (Pa)
temp = (data(:,6) + 459.67)/1.8; % Temperature (K)
test = regexp(direct, '[0-9]*', 'match');
u_scale = sqrt(2*R*temp(test)*(k*dp(test))/(k*(patm-sp(test)))); % Normalization velocity 
end