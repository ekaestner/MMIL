function ScanInfo = PING_anonymize(PING_dir)

ScanInfo.ErrMsg = '';

if ~exist(PING_dir), 
   ScanInfo.ErrMsg = sprintf('%s does not exist.',PING_dir); 
   return;
end

ScanInfo_file = sprintf('%s/ScanInfo.txt',PING_dir);
if exist(ScanInfo_file,'file')
   ScanInfo = PING_Proc_ScanInfoFile(ScanInfo_file);
   if ~isempty(ScanInfo.ErrMsg)
      return;
   end
   ScanInfo = PING_dir_dcmodify(PING_dir,ScanInfo);
%   [status,result]=system(sprintf('rm %s',ScanInfo_file));
%   if status, error(result); end
else
   cd(PING_dir);
   D=dir;
   for ii=1:length(D),
      if D(ii).name(1)=='.', continue; end
      if D(ii).isdir,
         subdirname = sprintf('%s/%s',PING_dir,D(ii).name);
         ScanInfo = PING_anonymize(subdirname);
      end
   end
end
