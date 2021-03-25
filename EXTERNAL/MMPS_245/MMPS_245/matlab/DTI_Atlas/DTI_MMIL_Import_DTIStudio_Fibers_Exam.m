function DTI_MMIL_Import_DTIStudio_Fibers_Exam(ContainerPath,FiberContainerPath,varargin);
%function DTI_MMIL_Import_DTIStudio_Fibers_Exam(ContainerPath,FiberContainerPath);
% 
% Purpose: Import the manual Fibers done in DTIStudio, then create count map for atlas creation
%
% Required Parameters:
%   ContainerPath: full path of DTIPROC Container
% 	FiberContainerPath: full path of Manual FiberContainer Path
%
% Optional Parameters:
%  'DTI_min_ndirs': minimum number of directions for ContainerPath
%    {default = 6}
%	 'DTI_min_bval': minimum b-value for ContainerPath
%    {default = 1000}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%	 'DTI_min_nb0': minimum number of b=0 images for ContainerPath
%    {default = 1}
%	 'DTI_snums': scan numbers to be included for ContainerPath
%    {default = []}
%	 'DTI_revflag': whether to include rev phase-encode scans for ContainerPath
%    {default = 0}
%	 'DTI_fiber_infix': DTI processed data infix for ContainerPath
%    {default = 'corr_regT1'}
% 
% Optional Parameters that determine manual fiber data
%  'Fiber_studytype': method to create fiber file path
%    0: use snums and infix
%    1: use min_bval and infix
%    {default = 0}
%  'Fiber_dirname': Name of manual fiber directory
%    {default = 'fibers'}
%  'Fiber_snums': scan numbers to be included for FiberContainerPath
%    {default = []}
%  'Fiber_infix': Fiber ContainerPath infix
%    {default = 'reg_mc_ecc_B0uw_gruw_iso'}
%  'Fiber_min_ndirs': minimum number of directions for FiberContainerPath
%    {default = 6}
%  'Fiber_min_bval': minimum b-value for FiberContainerPath
%    {default = 1000}
%  'Fiber_flex_flag': [0|1] DTI_flex included for FiberContainerPath
%    {default = 0}
%  'Fiber_min_nb0': minimum number of b=0 images for FiberContainerPath
%    {default = 1}
%	 'Fiber_DTIStudio_fpat': manual fibers filename pattern
%    {default = 'Fiber_(?<fnum>\d+)_path.dat'}
%	 'Fiber_revsliceflag': [0|1] reverse slice order to LPS  
%    {default = 1}
%	 'Fiber_permvec': permutation order in case the slices are not axial
%    {default = [1,2,3]}
%	 'Fiber_mask_datatype': datatype for mask output
%    {default = 'uint8'}
%	 'Fiber_count_datatype' datatype for count output
%    {default = 'float'}
%	 'Fiber_countflag': [0|1] output option as mask or count maps
%    {default = 1}
%
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Output Files:
%  Creates output files of mask and count maps of all the Fibers 
%
% Created:  02/17/11 by Vijay Venkatraman (Original code from Don Hagler)
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin,{...
	... % DTI ContainerPath parameters
	'DTI_min_ndirs',25,[],...
	'DTI_min_bval',1000,[],...
  'DTI_flex_flag',false,[false true],...
	'DTI_min_nb0',1,[],...
	'DTI_snums',[],[],...
	'DTI_revflag',2,[0 1 2],...
	'DTI_fiber_infix','corr_regT1',[],...
	... % Fiber Container
  'Fiber_studytype',0,[0 1]...
  'Fiber_dirname','fibers',[],...  
  'Fiber_snums',[],[],...
  'Fiber_infix','reg_mc_ecc_B0uw_gruw_iso',[],...
  'Fiber_min_bval',1000,[],...
  'Fiber_flex_flag',false,[false true],...
  'Fiber_min_nb0',1,[],...
  'Fiber_min_ndirs',6,[],...
	'Fiber_DTIStudio_fpat','Fiber_(?<fnum>\d+)_path.dat',[],...
	'Fiber_revsliceflag',1,[],...
	'Fiber_permvec',[1,2,3],[],...
	'Fiber_mask_datatype','uint8',[],...
	'Fiber_count_datatype','float',[],...
  'Fiber_countflag',1,[0 1],...
  ...
	'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ContainerInfo, get DTI scan info, determine valid scans, reference scans
%tags = {'snums','revflag','min_nb0','min_ndirs','min_bval'};
%args = mmil_parms2args(parms,tags);
%[ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
%if errcode ~= 0, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that DT calculations exist for ContainerPath
[DT_fstem,DT_snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,...
  'snums',parms.DTI_snums,'infix',parms.DTI_fiber_infix,...
  'min_ndirs',parms.DTI_min_ndirs,'min_bval',parms.DTI_min_bval,...
  'flex_flag',parms.DTI_flex_flag,'min_nb0',parms.DTI_min_nb0);
fname_dtmeas_in = [DT_fstem '_meas.mat'];
if ~exist(fname_dtmeas_in,'file')
	fprintf('%s: WARNING: %s not found, run MMIL_Process_DTI first\n',...
		mfilename,fname_dtmeas_in);
  return;
end;
% get dimensions of FA image of ContainerPath
fname_fa_in = sprintf('%s_FA.mgz',DT_fstem);
if ~exist(fname_fa_in,'file')
	fprintf('%s: WARNING: %s not found, run MMIL_Process_DTI first\n',...
		mfilename,fname_fa_in);
  return;
end;
[vol,M,mr_parms,volsz] = fs_load_mgh(fname_fa_in,[],[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FiberPath 
FiberPath = sprintf('%s/DTIStudio_fiber_paths',FiberContainerPath);
if ~exist(FiberPath,'dir')
	fprintf('%s: ERROR: input fiber path %s not found\n',mfilename,FiberPath);
  return;
end;

FiberFAPath = sprintf('%s/DTcalc',FiberContainerPath);
if ~exist(FiberFAPath,'dir')
	fprintf('%s: ERROR: input fiber path %s not found\n',mfilename,FiberFAPath);
  return;
end;
 
% Partial Fix as there is no ContainerInfo.mat
if parms.Fiber_studytype
   Fiber_fstem_tmp = DTI_MMIL_Set_DT_fstem(ContainerPath,...
    'snums',parms.Fiber_snums,'infix',parms.Fiber_infix,..
    'min_ndirs',parms.Fiber_min_ndirs,...
    'min_bval',parms.Fiber_min_bval,...
    'flex_flag',parms.Fiber_flex_flag,...
    'min_nb0',parms.Fiber_min_nb0);
   Fiber_fstem= strrep(Fiber_fstem_tmp,'/proc/',['/' parms.Fiber_dirname '/']);
elseif ~parms.Fiber_studytype
  Fiber_fstem_tmp=  DTI_MMIL_Set_DT_fstem(ContainerPath,...
    'snums',parms.Fiber_snums,'min_bval',parms.Fiber_min_bval,...
    'flex_flag',parms.Fiber_flex_flag,...
    'infix',parms.Fiber_infix);
  Fiber_fstem = strrep(Fiber_fstem_tmp,'/proc/',['/' parms.Fiber_dirname '/']);
  Fiber_fstem = strrep(Fiber_fstem,['_minb' int2str(parms.Fiber_min_bval)],'');
else
  fprintf('%s: ERROR: Check your manual Fibers directory\n',...
    mfilename);
  return;
end
 
% get dimensions of FA image of FiberContainerPath
fname_Fiberfa_in = sprintf('%s_FA.mgz',Fiber_fstem);
if ~exist(fname_Fiberfa_in,'file')
	fprintf('%s: WARNING: %s not found, run MMIL_Process_DTI first\n',...
    mfilename,fname_Fiberfa_in);
	return;
end;
[volf,Mf,mr_parmsf,volszf] = fs_load_mgh(fname_Fiberfa_in,[],[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output directories
DTIStudio_OutPath = sprintf('%s/DTIStudio_fiber_masks',ContainerPath);
mmil_mkdir(DTIStudio_OutPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Register the FA of FiberContainerPath to ContainerPath
fprintf('registering FA in Fiber Container to FA in Processed Container...\n');
fname_b0_in = sprintf('%s_b0.mgz',DT_fstem);
fname_fa_mask= sprintf('%s_%s_brainmask.mgz',DT_fstem,parms.Fiber_dirname);
fname_fiber_mask= sprintf('%s_%s_brainmask.mgz',Fiber_fstem,parms.Fiber_dirname);
fname_fiberb0_in= sprintf('%s_b0.mgz',Fiber_fstem);
fname_fiberfa_masked= sprintf('%s_%s_FA_masked.mgz',Fiber_fstem,parms.Fiber_dirname);

if exist(fname_fiberb0_in, 'file')
  epi_brainmask(fname_b0_in,fname_fa_mask,'forceflag',parms.forceflag);
  fprintf('%s: using dilated Fiber B0 image mask to remove background noise of Fiber FA image \n',mfilename);
  epi_brainmask(fname_fiberb0_in,fname_fiber_mask,'forceflag',parms.forceflag);
  volmask_fib_dilated= mmil_dilate_mask(ctx_load_mgh(fname_fiber_mask));
  volmask_fib_dilated_mgh= ctx_ctx2mgh(volmask_fib_dilated);
  volmask_fib_dilated_mgh_th= (volmask_fib_dilated_mgh>0);
  volf_mask= volf .* volmask_fib_dilated_mgh_th;
  fs_save_mgh(volf_mask,fname_fiberfa_masked,Mf);
  M_vol_to_volf = mmil_reg(fname_fa_in,fname_fiberfa_masked,'fname_maskA',fname_fa_mask,...
                          'outdir',FiberContainerPath,'affine_flag',0,'forceflag',parms.forceflag); % Note: M stores volf to vol
  cmd= sprintf('mv %s %s',fname_fa_mask,FiberContainerPath);
  [status,result]=unix(cmd);
  cmd= sprintf('mv %s %s',fname_fiber_mask,FiberContainerPath);
  [status,result]=unix(cmd);
else
  fprintf('%s:Using the original fiber FA image \n', mfilename);
  epi_brainmask(fname_b0_in,fname_fa_mask,'forceflag',parms.forceflag);
  M_vol_to_volf = mmil_reg(fname_fa_in,fname_Fiberfa_in,'fname_maskA',fname_fa_mask,...
                 'outdir',FiberContainerPath,'affine_flag',0,'forceflag',parms.forceflag); % Note: M stores volf to vol
  cmd= sprintf('mv %s %s',fname_fiber_mask,FiberContainerPath);
  [status,result]=unix(cmd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply the M_vol_to_volf to Fiber FA, can be used for QC
fname_fib2sub_resam= sprintf('%s/%s_fiber2subj_FA.mgz',DTIStudio_OutPath,parms.Fiber_dirname);
if exist (fname_fiberfa_masked,'file')
  [volfa_fib2sub_resam Mfib2sub_resam]= ctx_ctx2mgh(vol_resample_pad(...
              ctx_mgh2ctx(volf_mask,Mf),ctx_mgh2ctx(vol,M),M_vol_to_volf,2)); % cubic interpm
else
  [volfa_fib2sub_resam Mfib2sub_resam]= ctx_ctx2mgh(vol_resample_pad(...
              ctx_mgh2ctx(volf,Mf),ctx_mgh2ctx(vol,M),M_vol_to_volf,2)); % cubic interpm
end
fs_save_mgh(volfa_fib2sub_resam,fname_fib2sub_resam,Mfib2sub_resam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flist = dir(sprintf('%s/F*',FiberPath));
fprintf('%s: resampling Fibers from %s to %s...\n',...
   mfilename,FiberPath,DTIStudio_OutPath);
for f=1:length(flist)
	fname = char(flist(f).name);
	ftype_DTIStudio = regexp(fname,parms.Fiber_DTIStudio_fpat,'names');
	if ~isempty(ftype_DTIStudio)
		[fpath,fstem,fext] = fileparts(fname);
		if ~isempty(fext) & ~strcmp(fext,'.dat'), continue; end;
		fnum = str2double(ftype_DTIStudio.fnum);
		fname_in = sprintf('%s/%s',FiberPath,fname);
    fprintf('processing %s .. \n',fname);
    if ~isempty(parms.Fiber_countflag)
      if parms.Fiber_countflag, outfix= 'count';else, outfix= 'mask';end;
      fname_dat = sprintf('%s/fiber_%03d_%s.dat',DTIStudio_OutPath,fnum,outfix);
      fname_mat = sprintf('%s/fiber_%03d_%s.mat',DTIStudio_OutPath,fnum,outfix);
      if parms.Fiber_countflag | parms.forceflag 
        if ~exist(fname_dat,'file') | ~exist(fname_mat,'file')
          % convert Fiber path to Fiber count per voxel
          dti_fiber_path_to_mask(fname_in,fname_dat,1);
          if ~exist(fname_dat,'file')
            fprintf('%s: ERROR: failed to create Fiber count file %s\n',...
              mfilename,fname_dat);
            return;
          end;
          vol_Fiber = dti_load_dat(fname_dat,volszf,parms.Fiber_permvec,...
            parms.Fiber_revsliceflag,parms.Fiber_count_datatype);
				  [vol_Fiber_resample Mresam]= ctx_ctx2mgh(vol_resample_pad(...
              ctx_mgh2ctx(vol_Fiber,Mf),ctx_mgh2ctx(vol,M),M_vol_to_volf,0)); % NN interpm
          mmil_save_sparse(vol_Fiber_resample,fname_mat,M);
        end;
      elseif ~parms.Fiber_countflag | parms.forceflag
        if ~exist(fname_dat,'file')
          dti_fiber_path_to_mask(fname_in,fname_dat);
          if ~exist(fname_dat,'file')
            fprintf('%s: ERROR: failed to create Fiber mask file %s\n',...
               mfilename,fname_dat);
            return;
          end;
			  	vol_Fiber = dti_load_dat(fname_dat,volszf,parms.Fiber_permvec,...
            parms.Fiber_revsliceflag,parms.Fiber_mask_datatype);
          vol_Fiber(vol_Fiber>0)=1;
				  [vol_Fiber_resample Mresam]= ctx_ctx2mgh(vol_resample_pad(...
              ctx_mgh2ctx(vol_Fiber,Mf),ctx_mgh2ctx(vol,M),M_vol_to_volf,0)); % NN interpm
          mmil_save_sparse(vol_Fiber_resample,fname_mgh,M);
        end;
      end;
    end;
	end;
end;
