function DTI_MMIL_Oblique_Test(RootDir,varargin)
% function DTI_MMIL_Oblique_Test(RootDir,varargin)
%
% Required Input Parameters:
%   RootDir: full path of directory containing DTIPROC containers
%
% Optional Input Parameters:
%   'fname_out' : output filename (.csv)
%     {default = 'DTI_ObliqueList.csv'}
%   'fnamestem': file name stem used in DTI data
%     {default = 'DTI'}
%   'ext': file extension of processed DTI data
%     {default = 'mgz'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   obl_flag= Indicates whether the data is oblique
%   obl_angle= Indicates the oblique angle
%   Output will be saved as an excel sheet in the PROC Container
%
% Created:  08/03/11 by Vijay Venkatraman
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out','DTI_ObliqueList.csv',[],...
  'fnamestem','DTI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
});

col_labels = {...
  'VisitID','StudyDate','ScanNumber',...
  'Oblique_Flag','Oblique_AngleX',...
  'Oblique_AngleY','Oblique_AngleZ'};
separator= ',';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(RootDir,'dir') 
  error('RootDir %s not found',RootDir);
end;

if exist(parms.fname_out,'file') && ~parms.forceflag
  return;
end;

dirlist = dir(sprintf('%s/*',RootDir));

fid = fopen(parms.fname_out,'w');
for j=1:length(col_labels)
  if j>1, fprintf(fid,'%s',separator); end
    fprintf(fid,'"%s"',col_labels{j});
  end;
fprintf(fid,'\n');

fprintf('%s: checking whether DTI scans are oblique\n',mfilename);
for i=1:length(dirlist)
   ContainerDir = char(dirlist(i).name);
   regpat = '^DTIPROC_(?<VisitID>\w+)_(?<StudyDate>\d{8})\..+';
   n = regexp(ContainerDir,regpat,'names','once');
   if isempty(n) || ~dirlist(i).isdir
    continue;
   end;
   fprintf('%s: checking %s...\n',mfilename,ContainerDir);
   VisitID = n.VisitID;
   StudyDate = n.StudyDate;
   ContainerPath = [RootDir '/' ContainerDir]; 
   [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath);
   if errcode ~= 0, return; end;
   if ~isempty(SessInfo)
    if ~isempty(SessInfo.snums_DT)
      count = SessInfo.snums_DT;
    else
      count = SessInfo.snums_valid; 
    end;
   else
    fprintf('WARNING: SessInfo is empty for %s \n',ContainerDir);
    continue;
   end;
   for s = count
     fstem = sprintf('%s/%s%d',ContainerPath,parms.fnamestem,s);
     fname = sprintf('%s%s',fstem,parms.ext);
     if ~exist(fname,'file')
      continue;
     end;
     [obl_flag,obl_angle]=mmil_check_oblique(fname);
     fprintf(fid,'"%s"',VisitID); fprintf(fid,'%s',separator);
     fprintf(fid,'"%s"',StudyDate); fprintf(fid,'%s',separator);
     fprintf(fid,'"%d"',s); fprintf(fid,'%s',separator);
     fprintf(fid,'"%f"',obl_flag);  fprintf(fid,'%s',separator);
     fprintf(fid,'"%f"',obl_angle(1)); fprintf(fid,'%s',separator);
     fprintf(fid,'"%f"',obl_angle(2)); fprintf(fid,'%s',separator);
     fprintf(fid,'"%f"',obl_angle(3)); fprintf(fid,'\n');
   end;
end;   
fclose(fid);

return;
