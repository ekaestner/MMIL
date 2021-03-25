function rc_calc_fieldsign(subj,varargin)
%function rc_calc_fieldsign(subj,varargin)
%
% Usage:
%  rc_calc_fieldsign(subj,'key1', value1,...);
%
% Required parameters:
%   subj:  a string specifying the subject name (freesurfer recon dir)
%
% Optional parameters:
%  'polstem' - file stem for polar angle retinotopy data
%     e.g. for files polret_r-lh.mgh and polret_i-lh.mgh, polstem = "polret"
%     {default = 'polret'}
%  'eccstem' - file stem for eccentricity retinotopy data
%     e.g. for files eccret_r-lh.mgh and eccret_i-lh.mgh, eccstem = "eccret"
%     {default = 'eccret'}
%  'outstem'  - output file stem (omit extension, hemi)
%    {default = 'fieldsign'}
%  'poldir' - input directory for polar angle data
%     {default = pwd}
%  'eccdir' - input directory for eccentricity data
%     {default = pwd}
%  'outdir' - output directory for fieldsign calculations
%     {default = pwd}
%  'hemilist' - cell array containing list of hemispheres to process
%    {default = {'lh','rh'}}
%  'intype' - input file type ('mgh' or 'w')
%    {default: 'mgh'}
%  'smooth' - smoothing steps on surface before fieldsign calculation
%    {default: 50}
%  'revfs_flag' - [0|1] whether to reverse fieldsign
%    (switches colors on display)
%    {default: 1}
%  'surf' - surface to load
%    {default: white}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'forceflag' - [0|1] whether to overwrite existing output files
%    {default: 0}
%
% Created:  10/09/08 by Don Hagler
% Last Mod: 02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms( varargin, {...
  'polstem','polret',[],...
  'eccstem','eccret',[],...
  'outstem','fieldsign',[],...
  'poldir',pwd,[],...
  'eccdir',pwd,[],...
  'outdir',pwd,[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'intype','mgh',{'mgh','w'},...
  'smooth',50,[0,1000],...
  'revfs_flag',true,[false true],...
  'surf','white',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
... % hidden parms
  'infixlist',{'_r','_i'},[],...
  'cmdname','fieldsign',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if ~isempty(parms.subjdir)
  setenv('SUBJECTS_DIR',parms.subjdir)
else
  error('FreeSurfer SUBJECTS_DIR environment variable not set');
end;

if ~iscell(parms.hemilist), parms.hemilist = {parms.hemilist}; end;

if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create output dir %s:\n%s',parms.outdir,msg);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  % convert mgh to w
  for j=1:length(parms.infixlist)
    infix = parms.infixlist{j};
    if strcmp(parms.intype,'mgh') % must convert to w
      % convert to w for pol
      wfilestem = sprintf('%s/%s%s',...
        parms.outdir,parms.polstem,infix);
      wfile = [wfilestem '-' hemi '.w'];
      if ~exist(wfile,'file') || parms.forceflag
        mghfile = sprintf('%s/%s%s-%s.mgh',...
          parms.poldir,parms.polstem,infix,hemi);
        if ~exist(mghfile,'file'), error('file %s not found',mghfile); end;
        fs_mgh2w(mghfile,wfilestem,hemi,1);
      end;
      if ~exist(wfile,'file'), error('file %s not found',wfile); end;
      % repeat for ecc
      wfilestem = sprintf('%s/%s%s',...
        parms.outdir,parms.eccstem,infix);
      wfile = [wfilestem '-' hemi '.w'];
      if ~exist(wfile,'file') || parms.forceflag
        mghfile = sprintf('%s/%s%s-%s.mgh',...
          parms.eccdir,parms.eccstem,infix,hemi);
        if ~exist(mghfile,'file'), error('file %s not found',mghfile); end;
        fs_mgh2w(mghfile,wfilestem,hemi,1);
      end;
      if ~exist(wfile,'file'), error('file %s not found',wfile); end;
    else
      % check pol data exists
      wfilestem = sprintf('%s/%s%s',...
        parms.poldir,parms.polstem,infix);
      wfile = [wfilestem '-' hemi '.w'];
      if ~exist(wfile,'file'), error('file %s not found',wfile); end;
      % check ecc data exists    
      wfilestem = sprintf('%s/%s%s',...
        parms.eccdir,parms.eccstem,infix);
      wfile = [wfilestem '-' hemi '.w'];
      if ~exist(wfile,'file'), error('file %s not found',wfile); end;
    end;
  end;
  if strcmp(parms.intype,'mgh')
    parms.poldir = parms.outdir;
    parms.eccdir = parms.outdir;
  end;

  cmd = parms.cmdname;
  cmd = sprintf('%s -name %s',...
    cmd,subj);
  cmd = sprintf('%s -polstem %s -eccstem %s -outstem %s',...
    cmd,parms.polstem,parms.eccstem,parms.outstem);
  cmd = sprintf('%s -poldir %s -eccdir %s -outdir %s',...
    cmd,parms.poldir,parms.eccdir,parms.outdir);
  cmd = sprintf('%s -hemi %s -surf %s -smooth %d',...
    cmd,hemi,parms.surf,parms.smooth);
  if parms.revfs_flag
    cmd = [cmd ' -revfs'];
  end;    

  disp(cmd);
  [status,result] = unix(cmd);
  if status
    error('fieldsign failed:\n%s',result);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
