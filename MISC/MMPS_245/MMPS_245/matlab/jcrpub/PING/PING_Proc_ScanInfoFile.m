function [ScanInfo] = PING_Proc_ScanInfoFile(scan_info_file)

ScanInfo = new_ScanInfo(scan_info_file);   % subfunction

if ~exist(scan_info_file,'file')
   ScanInfo.ErrMsg = sprintf('%s file not found.',scan_info_file);
   return;
end

[fid,errmsg]=fopen(scan_info_file,'r');
if fid==-1, 
   ScanInfo.ErrMsg = errmsg;
   return;
end

delims = ': ';
while(1)
   tline = fgetl(fid);
   if ~ischar(tline), break, end

   [T,R] = strtok(tline,delims);
   switch T
     case {'SubjectID','VisitID','AttemptID','SiteID','ScheduledScanDate','DOB','UploadStatus','Sex'}
       [K,R] = strtok(R,delims);
       if isempty(K)
          ScanInfo.ErrMsg = sprintf('Missing %s value.',T);
          break;
       end
       ScanInfo = setfield(ScanInfo,T,K);
   end 
end
fclose(fid);

if strncmp(ScanInfo.UploadStatus,'Successful',10)
   ScanInfo.UploadStatus = 1;
else
   ScanInfo.UploadStatus = 0;
   ScanInfo.ErrMsg = sprintf('Upload process not successful.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaninfo = new_ScanInfo(filename)

scaninfo.FileName = filename;
scaninfo.SiteID = '';
scaninfo.SubjectID = '';
scaninfo.VisitID = '';
scaninfo.AttemptID = '';
scaninfo.ScheduledScanDate = '';
scaninfo.DOB = '';
scaninfo.Sex = '';
scaninfo.UploadStatus = '';
scaninfo.ErrMsg = '';
scaninfo.AnonOK = 0;
scaninfo.AnonFileCount = 0;
