function Re0 = Re0_gen_airfoil(direct)
airfoil_path = [direct, filesep 'Other' filesep 'A2_CpPHF_NS_20100728b.xls'];
[data, ~] = xlsread(airfoil_path, 'DataSheet');
test = regexp(direct, '[0-9]*', 'match');
if isempty(test)
    test = 2;
else
    test = str2double(test) + 2;
end
Re0 = data(:,2)*1000;
Re0 = Re0(test);
end