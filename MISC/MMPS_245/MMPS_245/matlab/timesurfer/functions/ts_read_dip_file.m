function dipole_info=ts_read_dip_file(dip_file)
% ts_read_dip_file.m:	read dipole file for calculating forward solution
%
% dipole_info=ts_read_dip_file(dip_file)
%
% read freesurfer dip file for the location and orientation of the dipoles
% returns dipole_info containing x,y,z location and nx,ny,nz normal vector
%
% created:    05/07/05 by Mingxiong Huang as read_dipdec.m
% early mod:  11/27/05 by Don Hagler as read_dip_file.m
% last mod:   08/01/06 by Don Hagler
%

dipole_info = [];

% make sure file exists
[status,message,messageid]=fileattrib(dip_file);
if (status==0)
  fprintf('%s: file %s not found... quitting\n',mfilename,dip_file);
  return;
end

% read dip file to get the location of the dipole
fprintf('ts_read_dip_file: reading file %s...\n',dip_file);
fp=fopen(dip_file,'r','ieee-be.l64');
[ch]=fread(fp,1,'uchar');
[n_dipole]=fread(fp,1,'int32');
dipole_info=fread(fp,[6,n_dipole],'float');
fclose(fp);

fprintf('%s: %d dipoles\n',mfilename,n_dipole);

