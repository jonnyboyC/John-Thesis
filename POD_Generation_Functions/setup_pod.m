problem_base1.num_images = 2000;
problem_base1.load_raw = false;
problem_base1.save_pod = true;
problem_base1.direct = 'D:\thesis\PIVData\Baseline';
problem_base1.save_figures = {'fig', 'jpg'};

problem_base2.num_images = 2000;
problem_base2.load_raw = false;
problem_base2.save_pod = true;
problem_base2.direct = 'D:\shear layer\PIVData\Baseline1';
problem_base2.save_figures = {'fig', 'jpg'};

problem_force1.num_images = 2000;
problem_force1.load_raw = false;
problem_force1.save_pod = true;
problem_force1.direct = 'D:\thesis\PIVData\Forced';
problem_force1.save_figures = {'fig', 'jpg'};

problem_force2.num_images = 2000;
problem_force2.load_raw = false;
problem_force2.save_pod = true;
problem_force2.direct = 'D:\shear layer\PIVData\Forced1';
problem_force2.save_figures = {'fig', 'jpg'};

problem_airfoil.num_images = 2000;
problem_airfoil.load_raw = false;
problem_airfoil.save_pod = true;
problem_airfoil.save_figures = {'fig', 'jpg'};
problem_airfoil.l_scale = 0.2032;
problem_airfoil.u_scale_gen = @u_scale_gen_airfoil;
problem_airfoil.flip_x = true;
problem_airfoil.flip_y = false;
