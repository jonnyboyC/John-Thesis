function [data,szones] = readgrid_old(fname,np,nps,data)
%data = READGRID(fname,[data])
%Reading grid file into the XG ,YG and ZG fields of DATA

%Read grid file FNAME to vector GRID

%fname = 'gri.dat';
fid = fopen(fname);
grid = fscanf(fid,'%g');%'double'

%JMAX, KMAX - number of points in each direction
%SZONE - matrix of dimentions of the zones
nzones = grid(1);
data.nzones = nzones;
jmax = grid(2);
kmax = grid(3);
szones = reshape(grid(2:2*nzones+1),2,nzones);

%Separate data part in GRID array and divide into 
%three cell arrays, DATA.XG, DATA.YG and DATA.ZG.
grid = grid(2*nzones+2:end);
[data.xg,data.yg,data.zg] = make3arr(grid,nzones,szones);
if nzones>1
for i=1:nzones
data.xg{i}=data.xg{i}(np,1:nps)*.3048;
data.zg{i}=data.zg{i}(np,1:nps)*.3048;
data.yg{i}=data.yg{i}(np,1:nps)*.3048;
end
else
data.xg{1}=data.xg{1}(np,:)*.3048;
data.zg{1}=data.zg{1}(np,:)*.3048;
data.yg{1}=data.yg{1}(np,:)*.3048;
end
%==========================================================%

function[x,y,z] = make3arr(grid,nzones,szones)
%Converting GRID vector into three cell structures
%XG, YG and ZG

%Initialize
x = cell(1, nzones);
y = cell(1, nzones);
z = cell(1, nzones);
offset = 0;

%Dividing Grid vector and  reshaping parts into matrices
for i = 1:nzones
   nx = szones(1,i);
   ny = szones(2,i);
   x{i} = reshape(grid([1:nx*ny]+offset),nx,ny);
   offset = offset+nx*ny;
   y{i} = reshape(grid([1:nx*ny]+offset),nx,ny);
   offset = offset+nx*ny;
   z{i} = reshape(grid([1:nx*ny]+offset),nx,ny);
   offset = offset+nx*ny;
   %if i == 1
   %   X=x{i};
   %   Y=y{i};
   %   Z=z{i};
   %else
   %   X=cat(2,X,x{i});
   %   Y=cat(2,Y,y{i});
   %   Z=cat(2,Z,z{i});
   %end
end

%==========================================================%
