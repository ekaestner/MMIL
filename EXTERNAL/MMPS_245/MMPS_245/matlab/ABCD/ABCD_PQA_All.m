function abcd_pqa_all(varargin)
%function abcd_pqa_all(varargin)
%
%
%  Created:                by Jose Teruel
%  Last Mod by: 09/28/2017 by Feng Xue
%%%%%%%%%%% Variables
if ~mmil_check_nargs(nargin,3), return; end;
% check input parameters
parms = check_input(varargin);

%%%%%%%%%%%%%%%% Get all folders %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% FOR INDIR 1

inputFoldersList = dir(parms.indir);
inputFoldersList = inputFoldersList([inputFoldersList.isdir] & ~strncmpi('.', {inputFoldersList.name}, 1)); 

%%%%%%%%%%%%% Observe non-study folders or already processed %%%%%

fname = fullfile(parms.outdir,'processedList.txt');
if exist(fname, 'file')
    finishedExam=fileread(fname);
else
    finishedExam= '';
end        
paramsfileID = fopen(sprintf('%s/batchdirs/%s_all/scriptlist.txt',getenv('HOME'),parms.ProjID),'w');
counter = 0;
for i=1:length(inputFoldersList)
    
    st_indir = fullfile(parms.indir,inputFoldersList(i).name);
    matchCheck=strfind(finishedExam, st_indir);
    
    if isempty(matchCheck)
      batchNameScriptlist = sprintf('job_PQA_%03d',i);
      batchName = sprintf('%s/batchdirs/%s_all/job_PQA_%03d.m',getenv('HOME'),parms.ProjID,i);
      st_outdir = fullfile(parms.outdir,inputFoldersList(i).name); 
      fprintf(paramsfileID,[batchNameScriptlist,'\n']);
      batchID = fopen(batchName,'w');
      fprintf(batchID,'abcd_pqa(''%s'', ''%s'')\nexit\n',st_indir,st_outdir);
      fclose(batchID);
      counter = counter + 1;
    end
end

fclose(paramsfileID);
fprintf('Number of jobs created: %d\n',counter);
return;
function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','/space/syn11/1/data/ABCD/ABCD_QA_proc',[],...
    'outdir','/space/syn11/1/data/ABCD/ABCD_QA_proc/Summary',[],...
    'ProjID','DAL_ABCD_PQA',[],...
  });
return;
