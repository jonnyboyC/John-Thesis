function data = readtec( fname,fields,data,np,nps)
%data = readtec( nzones, fields )
%Read one recplot file and include information into
%defined by user fields of DATA
%
%Written by Marina Khibnik 07/23/99

nzones = data.nzones;
%Read one tecplot file and return names and number 
%of variables and all the variable information
[nvar,var,zone] = read_tec(fname,nzones);

%NAME is vector of variables' names associated with user defined fields
for i = 1 : length(fields)
   name{i} = var{fields(i)};
end
data.name = name;

%Initialize fields of DATA
for j = 1 : length(name)
   cmd = ['data.' name{j} ' = cell(nzones,1);'];
   eval(cmd);
end

%Input information from tecplot file into the fields of DATA   
for z = 1: nzones
   for k = 1 : length(name)
      cmd = ['data.' name{k} '{z} = zone{z}.field(1:end,1:end,' num2str(fields(k)) ');'];
      eval(cmd);
      if nzones==1
         cmd = ['data.' name{k} '{z} = data.' name{k} '{z}(np,:);'];
      else
         cmd = ['data.' name{k} '{z} = data.' name{k} '{z}(np,1:nps);'];
      end
      eval(cmd);
   end
end

%Initialize X and Y fields of DATA
data.x = cell(1,nzones);
data.y = cell(1,nzones);
data.z = cell(1,nzones);

%Input information from tecplot file into the X and Y fields of DATA   
for z = 1 : nzones
%   data.x{z} = zone{z}.x(2:end-1,2:end-1);
%   data.y{z} = zone{z}.y(2:end-1,2:end-1);
   data.x{z} = zone{z}.x;
   data.y{z} = zone{z}.y;
   data.z{z} = zone{z}.z;
   if nzones==1
      data.x{z} = data.x{z}(np,:)*.3048;
      data.y{z} = data.y{z}(np,:)*.3048;
      data.z{z} = data.z{z}(np,:)*.3048;
   else
      data.x{z} = data.x{z}(np,1:nps)*.3048;
      data.y{z} = data.y{z}(np,1:nps)*.3048;
      data.z{z} = data.z{z}(np,1:nps)*.3048;
   end
end


%==================================================================

function [nvar,var,zone] = read_tec(fname,nzones)
% [NVAR,VAR,ZONE] = READ_TEC(FNAME,NZONES)
% read tecplot file and save data into Matlab structure.	
%
% FNAME is the name of the input tecplot file
% NVAR is the number of variables in the header of tecplot file. It is equal
%    to the number of fields plus two.
% VAR is the list of variables in the header of tecplot file. 
%    VAR has a structure of a cell array of strings, VAR=cell(1,NVAR).
% NZONES is the number of blocks (zones) in the file.
% ZONE contains data for each block. It is a (1 x NZONE) structure array 
% with fields:
%    nz 	- block number (or identifier)
%    nx		- number of grid points in x-direction
%    ny		- number of grid points in y-direction
%    x		- x-coordinates of grid points (size(x)=[nx ny])
%    y		- y-coordinates of grid points (size(y)=[nx ny])
%    field	- fields at grid points (size(fields)=[nx ny NVAR-2])
%
% Alexander I. Khibnik
% UTRC
% 17-Feb-98 
% Modified by Mike Dorobantu 5/3/99 to fit Van Slooten's format


nvar = 7;

fid = fopen(fname,'rt');
var = cell(1,nvar);

fgets(fid);
var{1} = 'x'; 
var{2} = 'y';
var{3} = 'z'; 
var{4} = 'r'; 
var{5} = 'u'; 
var{6} = 'v'; 
var{7} = 'w';

str = [];
for i=1:nvar
	str = [str ' %g'];
end

zone= cell(1,nzones);
iblk = 0;
while (1)
   clear line
   line = fgets(fid);
   if line == -1 | iblk==nzones; 
      break
   end 
   iblk = iblk+1;
	%nz=str2num(line(15:16));
   nx = str2num(line(15:17));
   ny = str2num(line(22:24));  
	x = zeros(nx,ny);
   y = zeros(nx,ny);
   z = zeros(nx,ny);
	field = zeros(nx,ny,nvar);
	for j =1:ny
      v=fscanf(fid, str, [nvar nx]);
      x(:,j)=v(1,:)';
      y(:,j)=v(2,:)';
      z(:,j)=v(3,:)';
      field(:,j,:)=v(:,:)';
   end
   zone{iblk} =struct('nx',nx,'ny',ny,'x',x,'y',y,'z',z,'field',field);
   line=fgets(fid);
end
fclose(fid);

%==================================================================


