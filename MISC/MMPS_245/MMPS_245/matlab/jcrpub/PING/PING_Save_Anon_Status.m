function [errmsg] = PING_Save_Anon_Status(anon_stat_file,ScanInfo)

errmsg = '';

DISP_ERROR_MSG = true;

if isempty(anon_stat_file) || ~ischar(anon_stat_file) || isempty(ScanInfo) || ~isstruct(ScanInfo),
   errmsg = sprintf('%s: at least one invalid input arg.',mfilename);
   if DISP_ERROR_MSG, fprintf(1,'%s\n',errmsg); end
   return;
end

[fid,errstr]=fopen(anon_stat_file,'w');
if fid==-1,
   errmsg = sprintf('%s: ERROR: %s',mfilename,errstr);
   if DISP_ERROR_MSG, fprintf(1,'%s\n',errmsg); end
   return;
end

fprintf(fid,'UploadStatus_File=%s\n',ScanInfo.FileName);
fprintf(fid,'SubectID=%s\n',ScanInfo.SubjectID);
fprintf(fid,'VisitID=%s\n',ScanInfo.VisitID);
fprintf(fid,'AttemptID=%s\n',ScanInfo.AttemptID);
fprintf(fid,'SiteID=%s\n',ScanInfo.SiteID);
fprintf(fid,'UploadOK=%d\n',ScanInfo.UploadStatus);
fprintf(fid,'AnonOK=%d\n',ScanInfo.AnonOK);
fprintf(fid,'AnonFileCount=%d\n',ScanInfo.AnonFileCount);
fprintf(fid,'Anon_ErrMsgLength=%d\n',length(ScanInfo.ErrMsg));
fprintf(fid,'Anon_ErrMsg=%s\n',ScanInfo.ErrMsg);
fprintf(fid,'EOF\n');
try fclose(fid);
catch
  errmsg = sprintf('%s: ERROR: problem writing file %s.',...
              mfilename,anon_stat_file);
  if DISP_ERROR_MSG, fprintf(1,'%s\n',errmsg); end
end
