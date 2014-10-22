%% This script simply changes the names of some variables in the workspace 
% so they comply with the naming convenctio of previous code development.
clc
clear all

load('D:\shear layer\PIVData\Data1\POD Data\POD.mat');
muz = {mean_u};
mvz = {mean_v};
phiutc = {pod_u1};
phivtc = {pod_v1};
vt = {vol_frac};
bi = bnd_idx;
sgtc2 = lambda2;
psitc = eig_func_norm;
x = {x};
y = {y};
sa = {dimensions};
t0 = [1 2];
Nsamples = 2000;
N = 42504;
