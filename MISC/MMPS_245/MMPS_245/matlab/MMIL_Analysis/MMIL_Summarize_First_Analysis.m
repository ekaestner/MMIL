function MMIL_Summarize_First_Analysis(ProjID,varargin)
% function MMIL_Summarize_First_Analysis(ProjID,[options])
%
% Required Input:
%   ProjID: project ID to for which to summarize analysis (e.g. 'REC_TEST')
%
% Optional Parameters
%  'StudyInfo': struct array of study information
%     (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields: 'proc' and 'fibers'
%     these specify the full paths of root directories containing data containers
%     {default = []}
%  'firstdir': where text file output from first processing
%     can be found placed relative to each ContainerPath
%     {default = 'firsthippo_nu'}
%  'firststem': output file stem from firsthippo.csh
%     {default = 'subj'}
%  'outdir': output directory
%     Created in user's home directory if not full path
%     if empty, will have name like: 'MetaData/ROI_Summaries' or
%                                    'MetaData/ProjID/ROI_Summaries'
%     {default = []}
%  'outstem': output file stem
%     {default = 'firsthippo'}
%  'qcflag': only include subjects with StudyInfo.QC=1
%     {default=1}
%  'forceflag': overwrite existing output
%     {default = 0}
%
%  Created:  11/02/11 by Don Hagler
%  Last Mod: 11/08/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'ProjID',ProjID,[],...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'firstdir','firsthippo_nu',[],...
  'firststem','subj',[],...
  'outdir',[],[],...
  'outstem','firsthippo_nu',[],...
  'qcflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'required_rootdirs',{'proc'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args= MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
N = length(StudyInfo);
if N==0
  fprintf('%s: ERROR: StudyInfo is empty',mfilename);
  return;
end;

if ~isempty(parms.ProjID)
  parms.outstem = [ProjID '_' parms.outstem];
end;

if isempty(parms.outdir)
  if ~isempty(parms.ProjID)
    parms.outdir = sprintf('%s/MetaData/%s/ROI_Summaries',...
      RootDirs.home,parms.ProjID);
  else
    parms.outdir = sprintf('%s/MetaData/ROI_Summaries',RootDirs.home);
  end;
end;
mmil_mkdir(parms.outdir);

fname_out = sprintf('%s/%s_volumes.csv',parms.outdir,parms.outstem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out,'file') || parms.forceflag
  fprintf('%s: summarizing FIRST analysis for %d studies...\n',...
    mfilename,N);
  fid=fopen(fname_out,'wt');
  wrote_header_flag = 0;
  nroi = 0;
  roinames = [];
  for i=1:N
    SubjID = StudyInfo(i).VisitID;
    VisitID = StudyInfo(i).VisitID;
    StudyDate = num2str(StudyInfo(i).StudyDate);
    ContainerPath = sprintf('%s/%s',RootDirs.proc,StudyInfo(i).proc);
    fname = sprintf('%s/%s/%s.volume',...
      ContainerPath,parms.firstdir,parms.firststem);
    if ~exist(fname,'file')
      fprintf('%s: WARNING: %s not found\n',mfilename,fname);
      continue; 
    end
    try
      first_stats = mmil_read_first_stats(fname);
    catch
      fprintf('%s: ERROR: failed to read first stats for VisitID %s:\n%s',...
        mfilename,VisitID,lasterr);
      continue;
    end;
    if ~wrote_header_flag
      nroi = length(first_stats);
      roinames = {first_stats.roiname};
      fprintf(fid,'"SubjID","VisitID","StudyDate"');
      for r=1:nroi
        fprintf(fid,',"FIRST_%s"',roinames{r});
      end;
      fprintf(fid,'\n');
      wrote_header_flag=1;
    end;
    tmp_roinames = {first_stats.roiname};
    if length(first_stats)~=nroi
      fprintf('%s: WARNING: VisitID %s has wrong number of ROIs in %s\n',...
        mfilename,VisitID,fname);
      continue;
    end;
    if length(tmp_roinames)~=length(intersect(tmp_roinames,roinames))
      fprintf('%s: WARNING: VisitID %s has mismatching ROIs in %s\n',...
        mfilename,VisitID,fname);
      continue;
    end;
    fprintf(fid,'"%s","%s","%s"',SubjID,VisitID,StudyDate);
    for r=1:nroi
      fprintf(fid,',%0.6f',first_stats(r).volume);
    end
    fprintf(fid,'\n');
  end;
  fclose(fid);
end;

