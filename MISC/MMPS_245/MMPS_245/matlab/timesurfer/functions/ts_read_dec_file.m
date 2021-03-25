function [dec_dipole,status]=ts_read_dec_file(dec_file)
% ts_read_dec_file.m:	read decimation file for calculating forward solution
%
% dec_dipole=ts_read_dec_file(dec_file)
%
% read binary freesurfer dec file for the dipoles (vertices) to include in the forward
%   solution
%
% created 05/07/05 by Mingxiong Huang as read_dipdec.m
% modified 11/27/05 by Don Hagler as read_dec_file.m
% last modified 08/01/06 by Don Hagler
%

if nargin<1
  help(mfilename);
  return;
end;

dec_dipole = [];
status = 0;
% make sure file exists
[status,message,messageid]=fileattrib(dec_file);
if (status==0)
  fprintf('%s: file %s not found... quitting\n',mfilename,dec_file);
  return;
end
% read dec file to get the indices of the decimated dipole.
%fprintf('ts_read_dec_file: reading file %s...\n',dec_file);
fp=fopen(dec_file,'r','ieee-be.l64');
[ch]=fread(fp,1,'uchar');
[n_dipole]=fread(fp,1,'int32');
[dec_dipole]=fread(fp,n_dipole,'uchar');
n_dec_dipole=length(find(dec_dipole));
fclose(fp);
%fprintf('%s: %d decimated dipoles\n',mfilename,n_dec_dipole);
status = 1;
return;
end
