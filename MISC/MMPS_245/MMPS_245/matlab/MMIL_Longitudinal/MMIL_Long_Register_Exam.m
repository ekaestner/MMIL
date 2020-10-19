function MMIL_Long_Register_Exam(dirA,dirB,varargin)
%function MMIL_Long_Register_Exam(dirA,dirB,[options])
%
% Purpose: Registers the T1 images for Longitudinal analysis
%
% Required Input:
%   dirA: full path of baseline directory
%   dirB: full path of followup directory
%
% Optional Parameters:
%  'nobiasflag':[0|1] register T1 using Quarc with no bias (forward and reverse)
%   {default: 1}
%  'forceflag': [0|1] if output files exist, delete them and then run setup
%   {default = 0}
%
% Created:  04/19/10 by Vijay Venkatraman
% Last Mod: 03/11/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin,{ ...
  'bindir',[],[],...
  'parmdir',[],[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'nobiasflag', true, [false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(parms.parmdir)
  parms.parmdir = [getenv('MMPS_PARMS') '/QUARC'];
end;

if isempty(parms.bindir)
  parms.bindir = [getenv('MMPS_DIR') '/bin'];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create script

[dirA_path,dirA_stem,dirA_ext] = fileparts(dirA);
outdir = [dirB '/nonlinreg_' dirA_stem,dirA_ext];

%if ~exist(outdir,'dir') || parms.forceflag
  fname_volA = [dirA '/images_resT1/T1' parms.ext];
  fname_volB = [dirB '/images_resT1/T1' parms.ext];
  fname_segA = [dirA '/images_resT1/aseg' parms.ext];
  fprintf('%s: registering T1 images...\n',mfilename);
  if parms.nobiasflag
    tic;  
    cmd = sprintf('quarcNoBias -b %s -f %s -a %s -bin %s -p %s -d %s',...
      fname_volA,fname_volB,fname_segA,parms.bindir,parms.parmdir,outdir);
    if parms.forceflag, cmd = [cmd '  -force']; end;       
    [status,result]= unix(cmd);
    if status
      error('%s: ERROR: cmd failed %s \n',mfilename,result);
    else
      fprintf('%s\n%s',cmd,result);    
    end;
   
    if parms.forceflag | ~exist([outdir '/f2b/T1_f2b_qc' parms.ext],'file')
      cmd = sprintf('applyTransforms -f %s -dx %s -dy %s -dz %s -m %s -t %s -cubic -o %s',...
        fname_volB,[outdir '/f2b/nonlinRegFine/dx' parms.ext],...
        [outdir '/f2b/nonlinRegFine/dy' parms.ext],...
        [outdir '/f2b/nonlinRegFine/dz' parms.ext],...
        [outdir '/f2b/affReg/affineRegMatrix.txt'],...
        fname_volA, [outdir '/f2b/T1_f2b_qc' parms.ext]);
           [status,result]= unix(cmd);
      if status
        error('%s:ERROR: cmd failed %s \n',mfilename,result);
      else
        fprintf('%s\n%s',cmd,result);    
      end;
    end;  
    seconds = toc;
    hours   = floor(seconds/3600);
    minutes = floor(((seconds-(3600*hours))/60));
    seconds = seconds-((3600*hours)+(60*minutes));
    fprintf('%s: Time elapsed: %02.0f:%02.0f:%02.0f.\n',...
      mfilename,hours,minutes,seconds);
  else
    tic;
    cmd= sprintf('quarc -b %s -f %s -a %s -bin %s -p %s -d %s',...
      fname_volA,fname_volB,fname_segA,parms.bindir,parms.parmdir,...
      [outdir '/f2b']);
    if parms.forceflag, cmd = [cmd '  -force']; end;         
    [status,result]= unix(cmd);
    if status
      error('%s:ERROR: cmd failed %s \n',mfilename,result);
    else
      fprintf('%s\n%s',cmd,result);    
    end;
    
    if parms.forceflag | ~exist([outdir '/f2b/T1_f2b_qc' parms.ext],'file')
      cmd= sprintf('applyTransforms -f %s -dx %s -dy %s -dz %s -m %s -t %s -cubic -o %s',...
        fname_volB,[outdir '/f2b/nonlinRegFine/dx' parms.ext],...
        [outdir '/f2b/nonlinRegFine/dy' parms.ext],...
        [outdir '/f2b/nonlinRegFine/dz' parms.ext],...
        [outdir '/f2b/affReg/affineRegMatrix.txt'],...
        fname_volA, [outdir '/f2b/T1_f2b_qc' parms.ext]);
      [status,result]= unix(cmd);
      if status
        error('%s:ERROR: cmd failed %s \n',mfilename,result);
      else
        fprintf('%s\n%s',cmd,result);    
      end;
    end;
    seconds = toc;
    hours   = floor(seconds/3600);
    minutes = floor(((seconds-(3600*hours))/60));
    seconds = seconds-((3600*hours)+(60*minutes));
    fprintf('%s: Time elapsed: %02.0f:%02.0f:%02.0f.\n',...
      mfilename,hours,minutes,seconds);
  end;
%end;

