function fs_write_wfile(fname, w, v);
%function fs_write_wfile(fname, w, v);
%
%  writes a vector into a binary 'w' file
%
%  fname - name of file to write to
%  w     - vector of values to be written
%  v     - vector of vertex indices for each value in w
%
% modified from freesurfer's fsfast toolbox fast_write_wfile.m
%
% created:  ?        by Don Hagler
% last mod: 08/20/07 by Don Hagler
%

if nargin<3
  help(mfilename);
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;

if fid==-1
  error('%s: failed to create w file %s\n',mfilename,fname);
end;

vnum = length(w) ;
if(vnum~=length(v))
  error('w and v must be same length');
end

fwrite(fid, 0, 'int16');
fs_fwrite3(fid, vnum);
for i=1:vnum
  k = v(i) - 1; % vertex numbers start at 0, but for matlab, index from 1
  fs_fwrite3(fid, k);
  wt = w(i);
  fwrite(fid, wt, 'float');
end

fclose(fid);

return

