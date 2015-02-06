setup_pod;
problem_airfoil.load_raw = true;
problem_base1.load_raw = true;
problem_force1.load_raw = true;
POD_Gen(problem_base1);
POD_Gen(problem_force1);
for i = 1:53
    problem_airfoil.direct = ['D:\thesis\PIVData\airfoil\Test_' num2str(i)];
    POD_Gen(problem_airfoil);
end

clear
setup_proj;
for i = 2:30
    problem_base1.num_modesG = i;
    Galerkin_Proj(problem_base1);
end

for i = 2:30
    problem_force1.num_modesG = i;
    Galerkin_Proj(problem_force1);
end


for i = 1:53
    problem_airfoil.direct = ['D:\thesis\PIVData\airfoil\Test_' num2str(i)];
    for j = 2:20
        problem_airfoil.num_modesG = j;
        Galerkin_Proj(problem_airfoil);
    end
end