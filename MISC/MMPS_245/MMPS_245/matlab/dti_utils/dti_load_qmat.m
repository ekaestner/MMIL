function qmat = dti_load_qmat(seq_type,nb0,ndirs,varargin)
%function qmat = dti_load_qmat(seq_type,nb0,ndirs,[options])
%
% Required Input:
%   seq_type: DTI sequence type
%    1: MGH Avanto(1.5T)/Allegra(3T) reversing PE polarity sequence
%    2: UCSD Symphony reversing PE polarity sequence
%    3: Other Siemens rorevro sequence
%    4: GE 12x
%    5: GE 14x
%    6: GE flex
%    11: NYU Siemens Allegra 3T
%    13: Other Siemens sequence
%    14: Philips
%  nb0: number b=0 scans
%  ndirs: number of diffusion gradient directions
%
% Optional Input:
%   'diffdirs': matrix of diffusion directions from dicoms
%     {default = []}
%   'tensor_fnum': tensor file number (for DTI_flex scans)
%     {default = 8}
%   'indir': input directory containing diffusion direction text files
%     if empty, will use $MMPS_parms/DTI
%     {default = []}
%
% Created:  01/26/09 by Don Hagler
% Last Mod: 02/02/17 by Don Hagler
%

qmat = [];

if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'diffdirs',[],[],...
  'tensor_fnum',8,[],...
  'indir',[],[],...
});

if isempty(parms.indir)
  MMPS_parms = getenv('MMPS_PARMS');
  if isempty(MMPS_parms), error('MMPS_PARMS environment variable not set'); end;
  parms.indir = [MMPS_parms '/DTI'];
end;
if ~exist(parms.indir,'dir'), error('directory %s not found',parms.indir); end;

switch seq_type
  case {1,2,3,11,13} % Siemens
    if isempty(parms.diffdirs)
      Siemens_DiffDirs_Prefix = [parms.indir '/Siemens_DiffDirs'];
      fname = sprintf('%s_%d.txt',Siemens_DiffDirs_Prefix,ndirs);
      if ~exist(fname,'file')
        error('diffusion direction file %s not found',fname);
      end;
      qmat = dti_read_dirs(fname,nb0);
      qmat(:,2) = -qmat(:,2);
    else
      qmat = parms.diffdirs;
    end; 
  case 4 % GE 12x
    if isempty(parms.diffdirs)
      fname = sprintf('%s/GE_tensor_12x.dat',parms.indir);
      if ~exist(fname,'file')
        error('diffusion direction file %s not found',fname);
      end;
      qmat = dti_read_GE_tensor_dat(ndirs,nb0,fname);
    else
      qmat = parms.diffdirs;
    end;
  case 5 % GE 14x
    if isempty(parms.diffdirs)    
      fname = sprintf('%s/GE_tensor_14x.dat',parms.indir);
      if ~exist(fname,'file')
        error('diffusion direction file %s not found',fname);
      end;
      qmat = dti_read_GE_tensor_dat(ndirs,nb0,fname);
    else
      qmat = parms.diffdirs;
    end;
  case 6 % GE DTI_flex
    fname = sprintf('%s/GE_flex_tensor%d.dat',...
      parms.indir,parms.tensor_fnum);
    if ~exist(fname,'file')
      error('diffusion direction file %s not found',fname);
    end;
    qmat = dti_read_GE_tensor_dat(ndirs,nb0,fname);
  case 14 % Philips
    if isempty(parms.diffdirs) || ndirs==32 
      fname = sprintf('%s/Philips_DiffDirs_32.txt',parms.indir);
      if ~exist(fname,'file')
        error('diffusion direction file %s not found',fname);
      end;
      qmat = dti_read_dirs(fname,nb0);
    else
      qmat = parms.diffdirs;
    end;
  otherwise
    error('unsupported DTI_Sequence_Type');
end;

