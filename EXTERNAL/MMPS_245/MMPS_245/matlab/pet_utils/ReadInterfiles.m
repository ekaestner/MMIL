function N = ReadInterfiles(fnames)

for i = 1:length(fnames)
  [N{i}.vol N{i}.M info] = interfileread_ADNI(fnames{i});
  N{i}.filetype = 'Interfile';
  N{i}.PatientID = info.PatientId;
  N{i}.StudyDateNum = datenum(sprintf('%s %s',info.StudyDateDdMmYryr,info.StudyTimeHhMmSs),'dd:mm:yyyy HH:MM:SS');
  N{i}.StudyDateStr = datestr(N{i}.StudyDateNum);
  N{i}.info = info; % Save rest of header, just in case...
end

