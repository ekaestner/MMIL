function fname_snr = mmil_fsl_calc_tSNR(fname,varargin)
%function fname_snr = mmil_fsl_calc_tSNR(fname,[options])
%
% Purpose:
% 	calculate temporal SNR using FSL tools
%
% Required Input:
% 	fname: full path of input file name (mgz or nii)
%
% Optional Parmaeters:
% 	'outdir': output directory
% 		{default = pwd}
%   'mcflag': [0|1] apply motion correction
%     {default = 1}
% 	'hp_sigma': high pass filter sigma (seconds)
% 		{default = 100}
% 	'TR': repetition time (seconds)
% 		{default = 1}
% 	'verbose': [0|1] display status messages
% 		{default = 0}
% 	'forceflag': [0|1] overwrite existing output
% 		{default = 0}
%
% Output:
% 	fname_snr: output file containing tSNR volume
%
% Created:  02/20/16 by Don Hagler
% Last Mod: 10/17/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_snr = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
	'outdir',[],[],...
  'mcflag',true,[false true],...
	'hp_sigma',100,[10,1000],...
	'TR',1,[0.01,100],...
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

if strcmp(parms.ext,'.nii')
  parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI; ');
elseif strcmp(parms.ext,'.nii.gz')
  parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI_GZ; ');
else
  parms.precmd = [];
end;

% calculate highpass sigma in volumes
parms.hp_sigma = parms.hp_sigma / parms.TR;

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

% motion correction
if parms.mcflag
  suffix = '_mcf';
  outstem = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix);
  fname_tmp = sprintf('%s%s',outstem,parms.ext);
  if ~exist(fname_tmp,'file') || parms.forceflag
	  if parms.verbose
		  fprintf('%s: performing motion correction...\n',mfilename);
		  tic;
	  end;
    cmd = sprintf('%smcflirt -in %s -out %s -refvol 0 -plots',...
								  parms.precmd,fname_in,outstem);
	  run_cmd(cmd);
	  if parms.verbose, toc; end;
  end;
  fname_in = fname_tmp;
else
  suffix = [];
end;

% high pass temporal filtering
suffix = [suffix '_hpf'];
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_tmp,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: applying high pass temporal filter...\n',mfilename);
		tic;
	end;
  cmd = sprintf('%sfslmaths %s -bptf %0.1f -1 %s',...
                parms.precmd,fname_in,parms.hp_sigma,fname_tmp);
	run_cmd(cmd);
	if parms.verbose, toc; end;
end;
fname_in = fname_tmp;

% calculate voxel-wise mean
suffix = '_mean';
fname_mean = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_mean,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating voxel-wise mean...\n',mfilename);
		tic;
	end;
  cmd = sprintf('%sfslmaths %s -Tmean %s',...
							  parms.precmd,fname_in,fname_mean);
	run_cmd(cmd);
	if parms.verbose, toc; end;
end;

% calculate voxel-wise standard deviation
suffix = '_std';
fname_std = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_std,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating voxel-wise standard deviation...\n',mfilename);
		tic;
	end;
  cmd = sprintf('%sfslmaths %s -Tstd %s',...
								parms.precmd,fname_in,fname_std);
	run_cmd(cmd);
	if parms.verbose, toc; end;
end;

% calculate temporal SNR
suffix = '_snr';
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_tmp,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating temporal SNR...\n',mfilename);
		tic;
	end;
  cmd = sprintf('%sfslmaths %s -div %s %s',...
							  parms.precmd,fname_mean,fname_std,fname_tmp);
	run_cmd(cmd);
	if parms.verbose, toc; end;
end;

% convert back to mgz
if ismember(fext,{'.mgh','.mgz'})
  fname_snr = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
	args = mmil_parms2args(parms,parms.convert_tags);
  suffix_list = {'_mean','_std','_snr'};
  for i=1:length(suffix_list)
    suffix = suffix_list{i};
    fname_in = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
    fname_out = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
    fs_mri_convert(fname_in,fname_out,args{:});
  end;
else
  fname_snr  = fname_tmp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_cmd(cmd)
	[s,r] = unix(cmd);
	if s, error('cmd %s failed:\n%s',cmd,r); end;
return;

