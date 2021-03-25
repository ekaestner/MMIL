function dti_convert_bvecs2dirs(fname_in,fname_out)
%function dti_convert_bvecs2dirs(fname_in,fname_out)
%
% Purpose: convert bvecs files from NYU Siemens DTI data to diffdirs
%
% Created:  05/26/09 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;


qmat = cell2mat(readtext(fname_in,' +'))';
fid = fopen(fname_out,'w');
for i=1:size(qmat,1)
  tmp = squeeze(qmat(i,:));
  if norm(tmp,1)
    fprintf(fid,'%f ',tmp);
    fprintf(fid,'\n');
  end;
end;
fclose(fid);

