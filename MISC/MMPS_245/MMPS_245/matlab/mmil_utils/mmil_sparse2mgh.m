function fname_mgh = mmil_sparse2mgh(fname_sparse,fname_mgh,forceflag)
%function fname_mgh = mmil_sparse2mgh(fname_sparse,[fname_mgh],[forceflag])
% 
% Purpose: convert sparse mat file to mgh format file
%
% Required Parameters:
%   fname_sparse: sparse input filename
%
% Optional Parameters:
%   fname_mgh: mgh output filename
%     if empty, will replace '.mat' extension with '.mgh'
%     {default = []}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_mgh: mgh output filename
%
% Created:  05/21/12 by Don Hagler
% Last Mod: 05/22/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('fname_mgh','var'), fname_mgh = []; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

[tpath,tstem,text] = fileparts(fname_sparse);

if ~strcmp(text,'.mat')
  error('input file does not have .mat extension');
end;

if isempty(fname_mgh)
  fname_mgh = regexprep(fname_sparse,'.mat','.mgh');
end;

if ~exist(fname_mgh,'file') || forceflag
  [vol,M] = mmil_load_sparse(fname_sparse);
  fs_save_mgh(vol,fname_mgh,M);
end;

return;

