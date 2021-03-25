function [L,geo]=load_tri_file(tri_file);
% function [L,geo]=load_tri_file(tri_file);
% It can handle two types of tri files:
% 
% Type 1: The tri file has 4 collumns [index x y z] for vartex location
%         and connection index also has an index for each triangle
% Type 2: The tri file has 4 collumns [x y z nx ny nz] for vartex location
%         and connection index has NO index for each triangle
%
% Output 
% L is n_ver by 6 [x y z nx ny nz] are the locations and normal direction 
%   of the vertices for Type1, no orientations for Typ2 2 
% geo is n_triangle by 3, the connection geometry matrix for triangles
%
% by M.X. Huang, PhD, March 2004.

fid=fopen(tri_file,'r');
line1=fgetl(fid);
n_ver=sscanf(line1,'%d');

line=fgetl(fid);
loc_temp=sscanf(line,'%f')';

if length(loc_temp)==4, % Type 1
    L=zeros(n_ver,3);
    L(1,:)=loc_temp(2:4); % the first line
    for i=2:n_ver, % the rest
        line=fgetl(fid);
        loc_temp=sscanf(line,'%f')';
        L(i,:)=loc_temp(2:4);
    end
elseif length(loc_temp)==6, % Type 2
    L=zeros(n_ver,6);
    L(1,:)=loc_temp; % the first line
    for i=2:n_ver, % the rest 
        line=fgetl(fid);
        L(i,:)=sscanf(line,'%f')';
    end
else
    error('Wrong tri file format');
end


% now the geo matrix
line_geo=fgetl(fid);
n_tri=sscanf(line_geo,'%d');

line=fgetl(fid);
geo_temp=sscanf(line,'%d')';

if length(geo_temp)==4, % Type 1
    geo=zeros(n_tri,3);
    geo(1,:)=geo_temp(2:4); % the first line
    for i=2:n_tri, % the rest
        line=fgetl(fid);
        geo_temp=sscanf(line,'%d')';
        geo(i,:)=geo_temp(2:4);
    end
elseif length(geo_temp)==3, % Type 2
    geo=zeros(n_tri,3);
    geo(1,:)=geo_temp; % the first line
    for i=2:n_tri, % the rest 
        line=fgetl(fid);
        geo(i,:)=sscanf(line,'%d')';
    end
else
    error('Wrong tri file format');
end


fclose(fid);