function REC_MMIL_Export_StudyInfo(StudyInfo,fname)
%function REC_MMIL_Export_StudyInfo(StudyInfo,fname)

% Created: 11/24/2008 AK
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
    fnames_bak=sprintf('%s/%s*',pathstr,name);
    [status,result] = unix(sprintf('cp -p %s %s',fnames_bak,dirname));
    if status
        error('failed to backup files %s to %s: %s',fnames_bak,dirname,result);
    end
    fprintf('%s: Finished backing up old files to: %s\n',mfilename,dirname);
end

% detect if it's a .csv or .mat file to be saved
[pathstr,name,ext] = fileparts(fname);
switch ext
  case '.csv'
    % loop thru fields in StudyInfo and output into spreadsheet
    headers = fieldnames(StudyInfo)';
    cols = length(headers);
    output = headers;
    for j=1:length(StudyInfo)
      data = struct2cell(StudyInfo(j));
      for k=1:cols
        output{j+1,k} = data{k};
      end
    end
    write_cellarray_csv(fname,output,[],0);
  case '.mat'
    save(fname,'StudyInfo');
end


