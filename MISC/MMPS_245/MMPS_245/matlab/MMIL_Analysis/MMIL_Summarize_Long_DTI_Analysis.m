function MMIL_Summarize_Long_DTI_Analysis(ProjID,varargin)
%function MMIL_Summarize_Long_DTI_Analysis(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%     (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%    'proc' and 'long_dti'
%     to specify the full paths of root directories containing data containers
%     {default = []}
%  'qcflag': only include subjects with StudyInfo.QC=1
%     {default = 1}
%
% Optional Parameters to select the DTI data
%  'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%  'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'nob0_flag': [0|1] whether to exclude b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_ndirs': minimum number of directions to be included in tensor fit 
%     {default = 6}
%  'min_nb0': minimum number of b=0 images required for tensor calculations (default=1)
%  'infix': infix of fibers used to warp to atlas (default= 'corr_regT1')
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%
% Optional Parameters:
%  'forceflag': create spreadsheets even if output files exist
%     {default = 0}
%
% Created:  05/10/11 by Vijay Venkatraman
% Last Mod: 07/07/15 by Don Hagler
%

% based on EPDTI_MMIL_PrePost_Summarize_Fibers, created 01/31/10 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters

FiberLegend = DTI_MMIL_Read_Fiber_Legend;
legend_fibernums = cell2mat({FiberLegend.FiberNumber});
legend_fibernames = {FiberLegend.FiberName};

parms = mmil_args2parms(varargin,{...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'qcflag',true, [false true],...
... % parameters for selecting the data
  'snums',[],[],...
  'min_bval',1000,[],...
  'flex_flag',false,[false true],...
  'nob0_flag',false,[false true],...
  'min_ndirs',6,[],...
  'min_nb0',1,[],...
  'infix','corr_regT1',[],...
  'revflag',0,[0,1,2],...
...  
  'forceflag',false,[false true],...
...
  'required_rootdirs',{'proc','long'},[],...
  'modality',{'MRI'},{'MRI','MEG','PET'},...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'resDTI_flag',false,[false true],...
  'dirprefix','LONG',[],...
  'outdir',[],[],...
  'measlist',{'FA','MD','TD','LD','b0N','T2w'},[],...
  'scalefacts',[1,1000,1000,1000,1],[],...
  'outfix',[],[],...
  'snums_flag',1,[0 1 2 3],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
...
  'fiber_outfix','fiber_roi_data',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get project specific info

tags_excl = {'StudyInfo','RootDirs','required_rootdirs','modality','qcflag', ...
  'snums','min_bval','flex_flag','nob0_flag','min_ndirs','min_nb0',...
  'infix','revflag','dirprefix'};
tags = setdiff(fieldnames(parms),tags_excl);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), return; end;
if isfield(StudyInfo(1),'Group'), SIflags.group = 1; else SIflags.group = 0;end;
if isfield(StudyInfo(1),'VisitCode'), SIflags.visitcode = 1; else SIflags.visitcode=0; end;
if isfield(StudyInfo(1),'Age'), SIflags.age = 1; else SIflags.age=0; end;
if isfield(StudyInfo(1),'StatsFlag'), SIflags.stats = 1;else SIflags.stats=0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify subjects with multiple visits

VisitNumbers = cell2mat({StudyInfo.VisitNumber});
ind = find(VisitNumbers>1);
SubjIDs = unique(({StudyInfo(ind).SubjID}));
ind = find(ismember({StudyInfo.SubjID},SubjIDs));
StudyInfo = StudyInfo(ind);
nsubs = length(SubjIDs);
if nsubs==0
  error('no subjects with followup scans found -- check StudyInfo VisitNumber');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.outfix)
  parms.outfix= ProjID;
end;

if parms.xcg_flag
  parms.outfix = [parms.outfix '_xcg'];
  parms.fiber_outfix = [parms.fiber_outfix '_xcg'];
end;
if parms.masksf_flag
  parms.outfix = [parms.outfix '_masksf'];
  parms.fiber_outfix = [parms.fiber_outfix '_masksf'];
end;

if parms.resDTI_flag
  dir_infix = 'resDTI';
else
  dir_infix = 'resT1';
end;

nmeas = length(parms.measlist);
if nmeas ~= length(parms.scalefacts)
  error('length of measlist (%d) does not match length of scalefacts (%d)',...
    nmeas,length(parms.scalefacts));
end;

parms.orig_measlist = parms.measlist;
parms.orig_scalefacts = parms.scalefacts;
parms.measlist = cell(1,3*nmeas);
parms.scalefacts = zeros(1,3*nmeas);
j = 1;
for i=1:length(parms.orig_measlist)
  meas = parms.orig_measlist{i};
  parms.measlist{j} = [meas '_vA'];
  parms.scalefacts(j) = parms.orig_scalefacts(i);
  j=j+1;
  parms.measlist{j} = [meas '_vB'];
  parms.scalefacts(j) = parms.orig_scalefacts(i);
  j=j+1;
  parms.measlist{j} = ['diff' meas];
  parms.scalefacts(j) = parms.orig_scalefacts(i);
  j=j+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.outdir)
  parms.outdir = sprintf('%s/MetaData/%s/ROI_Summaries',RootDirs.home,ProjID);
  mmil_mkdir(parms.outdir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: summarizing results for %d subjects...\n',mfilename,nsubs);
if ~isempty(parms.measlist)
  nrois = [];
  for m=1:length(parms.measlist)
    meas = parms.measlist{m};
    outfile = sprintf('%s/%s_%s.csv',parms.outdir,parms.outfix,meas);
    if ~exist(outfile,'file') || parms.forceflag
      fid=fopen(outfile,'wt');
      wrote_header_flag = 0;
      for i = 1:length(StudyInfo)
        SubjID = StudyInfo(i).SubjID;
        ind = find(strcmp(SubjID,{StudyInfo.SubjID}));
        tmpStudyInfo = StudyInfo(ind);
        ContainerPath = sprintf('%s/%s_%s',...
          RootDirs.long,parms.dirprefix,SubjID);
        vA = StudyInfo(i).VisitNumber;
        sessA = StudyInfo(i).SessID;
        dirA = sprintf('%s/visit%d',ContainerPath,vA);
        if ~exist(dirA,'dir')
          fprintf('%s: WARNING: dir %s not found',mfilename,dirA);
          continue;
        end;
        for j=1:length(tmpStudyInfo)
          vB = tmpStudyInfo(j).VisitNumber;
          if vB <= vA, continue, end;
          sessB = tmpStudyInfo(j).SessID;
          dirB = sprintf('%s/visit%d',ContainerPath,vB);
          if ~exist(dirB,'dir')
            fprintf('%s: WARNING: dir %s not found\n',mfilename,dirB);
            continue;
          end;
          indir = sprintf('%s/analysis_%s_visit%d',...
            dirB,dir_infix,vA);
          fname = [indir '/' meas '_' parms.fiber_outfix '.mat'];
          if ~exist(fname,'file')
            fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
            continue;
          end;
          fiber_data = [];
          load(fname);
          if isempty(fiber_data)
            fprintf('%s: WARNING: empty fiber_data in %s for %s\n',...
              mfilename,fname,SubjID);
            continue;
          end;
          tmp_nrois = length(fiber_data);
          if isempty(nrois)
            nrois = tmp_nrois;
          elseif tmp_nrois~=nrois
            % make sure nrois is same for each subject
            error('number of fiber ROIs for %s does not match first subject',...
              SubjID);
          end;
          if ~wrote_header_flag
            fprintf(fid,'"SubjectID","vA","vB"');
            % write additional fields if present in studyinfo
            if SIflags.group
              fprintf(fid,',"Group"');
            end;
            if SIflags.age
              fprintf(fid,',"Age"');
            end;
            if SIflags.stats
              fprintf(fid,',"StatsFlag"');
            end;
            for k=1:nrois
              if isfield(fiber_data(k),'fibernum')
                f = fiber_data(k).fibernum;
              else
                f = k;
              end;
              ind = find(f==legend_fibernums);
              if isempty(ind)
                fprintf(fid,',"Fiber_%d"',f);
              else
                fibername = strrep(legend_fibernames{ind},' ','-');
                fprintf(fid,',"DT_%s_%s_fibers"',meas,fibername);
              end;
            end;
            fprintf(fid,'\n');
            wrote_header_flag=1;
          end;
          fprintf(fid,'"%s",%d,%d',...
            SubjID,vA,vB);
          % write additional fields if present in studyinfo
          if SIflags.group
            fprintf(fid,',"%s"',StudyInfo(i).Group);
          end;
          if SIflags.age
            fprintf(fid,',%0.1f',StudyInfo(i).Age);
          end;
          if SIflags.stats
            fprintf(fid,',%d',StudyInfo(i).StatsFlag);
          end;
          scalefact = parms.scalefacts(m);
          for k=1:nrois
            if fiber_data(k).nvals_valid==0
              fprintf(fid,',NaN');
            else
              switch meas
                case 'vol'
                  % normalize by ICV
                  val = fiber_data(k).nvals/ICV;
                otherwise
                  val = fiber_data(k).avg;
              end;
              fprintf(fid,',%0.6f',...
                val*scalefact);
            end;
          end;
          fprintf(fid,'\n');
        end
      end;
      fclose(fid);
    end;
  end;
end;
