function REC_MMIL_Export_MRI_QC(StudyInfo,fname,dirlist_FS)
%function REC_MMIL_Export_MRI_QC(StudyInfo,fname,dirlist_FS)
%
% Created:  02/09/09 by Alain Koyama
% Rcnt Mod: 03/24/11 by Don Hagler
% Last Mod: 09/09/12 by Don Hagler
%

% back up existing output file
if exist(fname,'file')
  newdate=datestr(now, 'yyyymmdd');
  [pathstr, name, ext] = fileparts(fname);
  dirname=sprintf('%s/old/%s_Export/',pathstr,newdate);
  if ~exist(dirname,'dir')
    cmd = sprintf('mkdir -p %s',dirname);
    [status,result] = unix(cmd);
    if status
      error('failed to create backup directory %s: %s',dirname,result);
    end
  end
  fname_new=sprintf('%s/%s%s',dirname,name,ext);
  [status,result] = unix(sprintf('cp %s %s',fname,fname_new));
  if status
    error('failed to backup file %s to %s: %s',fname,fname_new,result);
  end
  fprintf('%s: Finished backing up old files to: %s\n',mfilename,dirname);
end

MRI_QC = mmil_readtext(fname,',','','"');

% Find relevant column in spreadsheet
headers = MRI_QC(1,:);
n = regexpi(headers,'^Subject ID$'); n = ~cellfun('isempty', n);
SubjID_col = find(n==1);

% extract Subject ID's and StudyDates into cells
DirNames = MRI_QC(2:end,SubjID_col);

tmpSubjIDs = regexp(DirNames,'(?<=FSURF_).+(?=_\d{8})','match');
SubjIDs = cell(1,length(tmpSubjIDs));
for i=1:length(SubjIDs), SubjIDs{i} = char(tmpSubjIDs{i}); end;
tmpStudyDates = regexp(DirNames,'(?<=_)\d{8}(?=\.)','match');
StudyDates = zeros(1,length(tmpStudyDates));
for i=1:length(StudyDates), StudyDates(i) = str2double(tmpStudyDates{i}); end;

for i=1:length(StudyInfo) % loop thru subjects in existing StudyInfo
  SubjID = StudyInfo(i).SubjID;
  StudyDate = StudyInfo(i).StudyDate;
  
  % find subject and visit in StudyInfo
  S_ind=find(strcmp(SubjIDs,SubjID) & StudyDates==StudyDate);
  if ~isempty(S_ind), continue; end; % subject already in qc file

  % if not found, add new entry in QC file if recon directory exists
  n = regexp({dirlist_FS.name},['FSURF_' SubjID '_' num2str(StudyDate) '.+'],'match');
  n(cellfun(@isempty,n)) = [];
  if isempty(n), continue; end; % no recon directory yet
  DirName = char(n{1});
  MRI_QC{end+1,SubjID_col} = DirName;

  % sort alphabetically
  headers = MRI_QC(1,:);
  MRI_QC_nohead = MRI_QC(2:end,:);
  [SubjIDs_sorted,idx] = sort(MRI_QC_nohead(:,1));
  MRI_QC = cat(1,headers,MRI_QC_nohead(idx,:));
end

write_cellarray_csv(fname,MRI_QC,[],0);
