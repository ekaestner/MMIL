function fname_brainmask = mmil_fsl_brainmask(fname,varargin)
%function fname_brainmask = mmil_fsl_brainmask(fname,[options])
%
% Purpose:
% 	create brain mask using FSL's brain extraction tool (bet)
%
% Required Input:
% 	fname: full path of input file name (mgz or nii)
%
% Optional Parmaeters:
% 	'outdir': output directory
% 		{default = pwd}
% 	'verbose': [0|1] display status messages
% 		{default = 0}
% 	'forceflag': [0|1] overwrite existing output
% 		{default = 0}
%
% Output:
% 	fname_brainmask: output file containing brain mask volume
%
% Created:  02/20/16 by Don Hagler
% Last Mod: 02/20/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_snr = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
	'outdir',[],[],...
	'verbose',false,[false true],...
	'forceflag',false,[false true],...
... % hidden
	'out_orient',[],[],...
	'ext','.nii',{'.nii','.nii.gz'},...
	'convert_tags',{'out_orient','verbose','forceflag'},[],...
});

if ~exist(fname,'file')
	error('file %s not found',fname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% convert to nii if necessary
[fpath,fstem,fext] = fileparts(fname);
switch fext
	case {'.mgh','.mgz'}
		fname_in = sprintf('%s/%s%s',parms.outdir,fstem,parms.ext);
		args = mmil_parms2args(parms,parms.convert_tags);
		fs_mri_convert(fname,fname_in,args{:});
	case {'.nii','.nii.gz'}
		fname_in = fname;
	otherwise
		error('invalid file extension: %s\n',mfilename,fext);
end;

% create brain mask
suffix = '_brainmask';
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_tmp,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: creating brainmask...\n',mfilename);
		tic;
	end;
  if strcmp(parms.ext,'.nii')
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI; ');
  elseif strcmp(parms.ext,'.nii.gz')
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI_GZ; ');
  else
    cmd = [];
  end;
  cmd = sprintf('%sbet %s %s',cmd,fname_in,fname_tmp);
%  cmd = sprintf('%sbet %s %s -m -B',cmd,fname_in,fname_tmp);
%% note: tried these options to correct for bias
%%       took longer, did not fix problem
	run_cmd(cmd);
	if parms.verbose, toc; end;
end;

% convert back to mgz
if ismember(fext,{'.mgh','.mgz'})
  fname_brainmask = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
	args = mmil_parms2args(parms,parms.convert_tags);
  fs_mri_convert(fname_tmp,fname_brainmask,args{:});
else
  fname_brainmask  = fname_tmp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_cmd(cmd)
	[s,r] = unix(cmd);
	if s, error('cmd %s failed:\n%s',cmd,r); end;
return;

