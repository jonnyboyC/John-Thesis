function Re0 = Re0_gen_airfoil(direct, ~, ~)
airfoil_path = [direct, '\Other\A2_CpPHF_NS_20100728b.xls'];
[data, ~] = xlsread(airfoil_path, 'DataSheet');
test = regexp(direct, '[0-9]*', 'match');
test = str2double(test);
Re0 = data(:,2)*1000;
Re0 = Re0(test);
end