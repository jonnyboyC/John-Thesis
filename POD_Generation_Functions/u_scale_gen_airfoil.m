function u_scale = u_scale_gen_airfoil(~, direct)
% U_SCALE_GEN_AIRFOIL pull velocity scale from file
%
%   u_scale = U_SCALE_GEN_AIRFOIL(~, direct)

airfoil_path = [direct filesep 'Other' filesep 'A2_CpPHF_NS_20100728b.xls'];
R       = 287;      % Gas constant
water2Pa = 249.17;  % inches water column to Pa
k       = .952;     % Tunnel constant
patm    = 98600;    % atmospheric pressure
data    = xlsread(airfoil_path, 'DataSheet'); % read values from xls sheet DataSheet
dp      = water2Pa*data(:,4)/k;     % Dynamic Pressure (Pa)
sp      = water2Pa*data(:,5);       % Static Pressure (Pa)
temp    = (data(:,6) + 459.67)/1.8; % Temperature (K)
test    = regexp(direct, '[0-9]*', 'match');
if isempty(test)
    test = 2;
else
    test = str2double(test) + 2;
end

% Calculate velocity scale based on xls spreadsheet values
u_scale = sqrt(2*R*temp(test)*(k*dp(test))/(k*(patm-sp(test))));  
end