problem_base1.num_images = 2000;
problem_base1.load_raw = false;
problem_base1.save_pod = true;
problem_base1.l_scale = 0.3048;
problem_base1.direct = 'D:\thesis\PIVData\shear layer\Baseline';
problem_base1.u_scale_gen = @u_scale_gen_shear;
problem_base1.save_figures = {'fig', 'jpg'};
problem_base1.flip = [true, false, true, false];

problem_base2.num_images = 2000;
problem_base2.load_raw = false;
problem_base2.save_pod = true;
problem_base2.l_scale = 0.3048;
problem_base2.direct = 'D:\shear layer\PIVData\shear layer\Baseline1';
problem_base2.u_scale_gen = @u_scale_gen_shear;
problem_base2.save_figures = {'fig', 'jpg'};
problem_base2.flip = [true, false, true, false];

problem_force1.num_images = 2000;
problem_force1.load_raw = false;
problem_force1.save_pod = true;
problem_force1.l_scale = 0.3048;
problem_force1.direct = 'D:\thesis\PIVData\shear layer\Forced';
problem_force1.u_scale_gen = @u_scale_gen_shear;
problem_force1.save_figures = {'fig', 'jpg'};
problem_force1.flip = [true, false, true, false];

problem_force2.num_images = 2000;
problem_force2.load_raw = false;
problem_force2.save_pod = true;
problem_force2.l_scale = 0.3048;
problem_force2.direct = 'D:\shear layer\PIVData\shear layer\Forced1';
problem_force2.u_scale_gen = @u_scale_gen_shear;
problem_force2.save_figures = {'fig', 'jpg'};
problem_force2.flip = [true, false, true, false];

problem_airfoil.num_images = 2000;
problem_airfoil.load_raw = false;
problem_airfoil.save_pod = true;
problem_airfoil.save_figures = {'fig', 'jpg'};
problem_airfoil.l_scale = 0.2032;
problem_airfoil.u_scale_gen = @u_scale_gen_airfoil;
problem_airfoil.flip = [false, false, true, false];

problem_cavity.num_images = 2000;
problem_cavity.load_raw = false;
problem_cavity.save_pod = true;
problem_cavity.save_figures = {'fig', 'jpg'};
problem_cavity.u_scale_gen = @u_scale_gen_shear;
problem_cavity.l_scale = 0.0127;

