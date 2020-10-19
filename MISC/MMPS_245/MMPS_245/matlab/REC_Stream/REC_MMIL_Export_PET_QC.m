function REC_MMIL_Export_PET_QC(StudyInfo,fname)
%function REC_MMIL_Export_PET_QC(StudyInfo,fname)

% Created 12/01/2008 AK
% Last modified: 06/10/2009 AK

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

PET_QC = mmil_readtext(fname,',','','"');

% Find relevant column in spreadsheet
headers = PET_QC(1,:);
n = regexpi(headers,'^SubjID$'); n = ~cellfun('isempty', n);
SubjID_col = find(n==1);
n = regexpi(headers,'^StudyDate_PET$'); n = ~cellfun('isempty', n);
StudyDate_col = find(n==1);

% extract Subject ID's and StudyDates into cells
SubjIDs = PET_QC(2:end,SubjID_col);
StudyDates = cell2mat(PET_QC(2:end,StudyDate_col));

% remove those subjects w/o a PET study date
SubjIDs_noPET = SubjIDs(cellfun('isempty',PET_QC(2:end,StudyDate_col)));
SubjIDs = SubjIDs(~cellfun('isempty',PET_QC(2:end,StudyDate_col)));

for i=1:length(StudyInfo) % loop thru subjects in existing QC file
  SubjID = StudyInfo(i).SubjID;
  StudyDate = StudyInfo(i).StudyDate_PET;
 
  % find subject and visit in StudyInfo
  if isempty(StudyDate) % just add blank extry w/ subjid if no PET
    S_ind=find(strcmp(SubjIDs_noPET,SubjID));
    if ~isempty(S_ind), continue; end; % if blank entry already exists in PET QC file, skip
    PET_QC{end+1,SubjID_col} = SubjID;
  else % if PET exists, locate entry in PET QC file
    S_ind=find(strcmp(SubjIDs,SubjID) & StudyDates==StudyDate);
    if ~isempty(S_ind), continue; end; % if already exists in PET QC file, skip
    % if not found, add new entry in QC file
    PET_QC{end+1,SubjID_col} = SubjID;
    PET_QC{end,StudyDate_col} = StudyDate;
  end; 

  % sort alphabetically
  headers = PET_QC(1,:);
  PET_QC_nohead = PET_QC(2:end,:);
  [SubjIDs_sorted,idx] = sort(PET_QC_nohead(:,1));
  PET_QC = cat(1,headers,PET_QC_nohead(idx,:));
end

write_cellarray_csv(fname,PET_QC,[],0);
