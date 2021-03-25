function T = fs_read_talxfm(fname)
%function T = fs_read_talxfm(fname)
%
% Purpose: read FreeSurfer talairach.xfm file
% 
% Input:
%   fname: FreeSurfer talairach.xfm file
%          found in <subject>/mri/transforms/talairach.xfm
%
% Output:
%   T: 3x4 matrix containing affine transformation
%     between FreeSurfer MRI volume and MNI template
%
% Note: Assuming surface RAS coordinates are in a vector,
%     v = [R A S 1]', and the talairach.xfm matrix is T,
%       then MNI Talairach coordinates are:
%       vMNI = T * v;
%
% Created:  11/17/04 by Darren Webber
% Last Mod: 11/03/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Revision: 1.1 $ $Date: 2004/11/17 21:03:27 $
% Licence:  GNU GPL, no implied or express warranties
% History:  11/2004 Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
T = [];

fid = fopen(fname,'r');
if isequal(fid,-1)
  error('failed to open file %s',fname);
else
  fprintf('%s: reading FreeSurfer talairach.xfm fname:\n%s\n',...
    mfilename,fname);

  % read lines until we get the string 'Linear_Transform', which precedes
  % the data transformation matrix
  success = 0;
  string2match = 'Linear_Transform';
  for i=1:20, % read up to 20 lines, no more
    tmp = fgetl(fid);
    if strmatch(string2match,tmp),
      % we have the right line, so don't read any more
      success = 1;
      break;
    end
  end

  if success
    % Read the transformation matrix (3x4).
    T = fscanf(fid,'%f',[4,3])';
    fclose(fid);
  else
    fclose(fid);
    error('failed to find ''Linear_Transform'' string in first 20 lines of xfm file.');
  end
end

return
